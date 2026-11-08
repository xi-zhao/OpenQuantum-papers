"""stories15M：模型里每一次"乘以权重"的运算，以及注意力里的 q·k 与 p·V，全部由量子干涉线路在 Qiskit Aer 上完成。

量子矩阵-向量乘（一条线路算一整个矩阵 M 乘向量 x）：
    |Ψ> = 1/sqrt(2R) Σ_j |j>_idx ⊗ (|0>_a |x̂>_data + |1>_a |m̂_j>_data)
    对辅助比特 a 做 H，测量 a 与行号 j：P(j,0) - P(j,1) = <x̂|m̂_j> / R
    y_j = |x| · |m_j| · cos_j
x̂ 与每一行 m̂_j 分处两个分支，乘积只在 H 门的干涉之后出现在测量统计里。

为了让"所有权重乘法"都进线路：
  - RMSNorm 的增益 g 并入后面的矩阵（W·diag(g)）
  - RoPE 旋转按位置并入 Wq、Wk（每个位置一组旋转后的行）
经典端只做：取词向量、向量长度、softmax、SiLU 与门控逐元素乘、残差相加、argmax、把计数换算成数值。

测量：Aer 给出线路末态在被测比特上的概率分布，再按 shots 数做多项分布抽样
（与 Aer 内部的 shots 采样是同一个过程，在 10^8 以上的 shots 量级快得多）；shots=None 时直接用精确概率。
"""
import math
import time
import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.library import SetStatevector

from model import load

SIM = AerSimulator(method="statevector")


def nbits(n):
    return max(1, math.ceil(math.log2(n)))


class QMatVec:
    """一组同时提交给 Aer 的量子矩阵-向量乘。"""

    def __init__(self, owner):
        self.o = owner
        self.jobs = []

    def add(self, x, M, Mnorm=None):
        self.jobs.append((x, M, np.linalg.norm(M, axis=1) if Mnorm is None else Mnorm))
        return len(self.jobs) - 1

    def run(self):
        circs, metas = [], []
        for x, M, Mn in self.jobs:
            R, n_in = M.shape
            ni, nd = nbits(R), nbits(n_in)
            D = 2 ** nd
            psi = np.zeros((2 ** ni, D, 2))
            xn = np.linalg.norm(x)
            psi[:R, :n_in, 0] = x / xn
            safe = np.where(Mn > 0, Mn, 1.0)
            psi[:R, :n_in, 1] = M / safe[:, None]
            psi = (psi.ravel() / math.sqrt(2 * R)).astype(complex)
            n = 1 + nd + ni
            qc = QuantumCircuit(n)
            qc.append(SetStatevector(psi), range(n))
            qc.h(0)
            qc.save_probabilities([0] + list(range(1 + nd, n)))
            circs.append(qc)
            metas.append((R, xn, Mn))
            self.o.stats["max_qubits"] = max(self.o.stats["max_qubits"], n)
        t = time.time()
        res = SIM.run(circs).result()
        self.o.stats["aer_seconds"] += time.time() - t
        outs = []
        for i, (R, xn, Mn) in enumerate(metas):
            p = res.data(i)["probabilities"].reshape(-1, 2)[:R]      # [j, a]
            if self.o.shots_per_row is None:
                n0, n1 = p[:, 0], p[:, 1]
            else:
                S = int(self.o.shots_per_row * R)
                cnt = self.o.rng.multinomial(S, (p / p.sum()).ravel()).reshape(R, 2)
                n0, n1 = cnt[:, 0].astype(float), cnt[:, 1].astype(float)
                self.o.stats["total_shots"] += S
            tot = n0 + n1
            cos = np.where(tot > 0, (n0 - n1) / np.maximum(tot, 1e-300), 0.0)
            outs.append(xn * Mn * cos)
        self.o.stats["circuits"] += len(circs)
        self.jobs = []
        return outs


class FullQuantumRunner:
    def __init__(self, shots_per_row=None, seed=0):
        self.cfg, W = load()
        c = self.cfg
        self.shots_per_row = shots_per_row
        self.rng = np.random.default_rng(seed)
        self.stats = dict(circuits=0, total_shots=0, max_qubits=0, aer_seconds=0.0)
        # RMSNorm 增益并入矩阵
        self.wq = W["wq"] * W["rms_att"][:, None, :]
        self.wk = W["wk"] * W["rms_att"][:, None, :]
        self.wv = W["wv"] * W["rms_att"][:, None, :]
        self.wo = W["wo"]
        self.w13 = np.concatenate([W["w1"], W["w3"]], axis=1) * W["rms_ffn"][:, None, :]
        self.w2 = W["w2"]
        self.cls = W["cls"] * W["rms_final"][None, :]
        self.cls_norm = np.linalg.norm(self.cls, axis=1)
        self.tok = W["tok"]
        i = np.arange(0, c.dim, 2)
        self.freq = 1.0 / (10000.0 ** ((i % c.hs) / c.hs))

    def rope_rows(self, Wm, pos):
        """把位置 pos 的 RoPE 旋转并入矩阵的行。"""
        cth, sth = np.cos(pos * self.freq)[:, None], np.sin(pos * self.freq)[:, None]
        out = np.empty_like(Wm)
        out[0::2] = cth * Wm[0::2] - sth * Wm[1::2]
        out[1::2] = sth * Wm[0::2] + cth * Wm[1::2]
        return out

    @staticmethod
    def rms_dir(x, eps=1e-5):
        return x / math.sqrt((x * x).mean() + eps)   # 只除以一个标量，增益已在矩阵里

    def generate(self, tokens, n_new, on_token=None):
        c = self.cfg
        L, H, hs, d = c.L, c.H, c.hs, c.dim
        Kc = [[] for _ in range(L)]; Vc = [[] for _ in range(L)]
        seq = list(tokens)
        pos = 0
        qmv = QMatVec(self)
        while True:
            x = self.tok[seq[pos]].copy()
            for l in range(L):
                # Q、K、V（含 RoPE 与 RMSNorm 增益）：一条线路
                M = np.concatenate([self.rope_rows(self.wq[l], pos), self.rope_rows(self.wk[l], pos), self.wv[l]])
                qmv.add(self.rms_dir(x), M)
                qkv, = qmv.run()
                q, k, v = qkv[:d], qkv[d:2 * d], qkv[2 * d:]
                Kc[l].append(k); Vc[l].append(v)
                K, V = np.array(Kc[l]), np.array(Vc[l])
                sl = [slice(h * hs, (h + 1) * hs) for h in range(H)]
                # 注意力打分 q·k：每个头一条线路
                for s_ in sl:
                    qmv.add(q[s_], K[:, s_])
                scores = [sc / math.sqrt(hs) for sc in qmv.run()]
                probs = []
                for sc in scores:
                    e = np.exp(sc - sc.max()); probs.append(e / e.sum())
                # 加权求和 p·V：每个头一条线路（行 = 头内 48 个维度）
                for h, s_ in enumerate(sl):
                    qmv.add(probs[h], V[:, s_].T)
                att = np.concatenate(qmv.run())
                # 输出投影 Wo
                qmv.add(att, self.wo[l])
                o, = qmv.run()
                x = x + o
                # 前馈：W1、W3 一条线路，W2 一条线路
                qmv.add(self.rms_dir(x), self.w13[l])
                h13, = qmv.run()
                h1, h3 = h13[:c.hidden], h13[c.hidden:]
                g = h1 / (1 + np.exp(-h1)) * h3
                qmv.add(g, self.w2[l])
                f, = qmv.run()
                x = x + f
            pos += 1
            if pos >= len(seq):
                if len(seq) - len(tokens) >= n_new:
                    return seq
                # LM head：32000 行，分 8 条线路
                xf = self.rms_dir(x)
                for a in range(0, c.vocab, 4000):
                    qmv.add(xf, self.cls[a:a + 4000], self.cls_norm[a:a + 4000])
                logits = np.concatenate(qmv.run())
                seq.append(int(logits.argmax()))
                if on_token:
                    on_token(seq)
