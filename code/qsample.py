"""任务 #23：把 QZ03 量子随机数池接入 stories15M 的温度/top-p 采样，与 MT19937、CSPRNG 对照。

设计（量子算力中心场景）：GPU/CPU 负责前向矩阵乘，QPU 只提供"解码采样所需的随机数"。
三种随机源都统一量化为 16-bit 均匀数 u∈[0,1)，公平对比：
  - QRNG ：QZ03 真机比特经 Toeplitz 提取后的池（有限，549600 bit）
  - MT19937：Python 标准库 random（梅森旋转，可预测——见 #24）
  - CSPRNG：操作系统密码学随机数（os.urandom）

实验 A（分布保真度）：对"little"后首个新 token 的 top-p 定分布 P，用三种源各抽 N 次，
  卡方拟合优度（经验计数 vs N·P）、消耗的 u 的 KS 均匀性、与 P 的总变差距离。
  → 验证 QRNG 是"无损即插即用"：都能复现目标分布。
实验 B（多样性）：各生成若干条续写，比较 distinct-n、唯一续写占比、token 频率相关性。
"""
import json
import os
import random
import numpy as np
from math import erf, sqrt
from scipy.stats import chisquare, kstest

import sys
sys.path.insert(0, "/tmp/claude-0/-home-claude/2b4b55be-e4b8-52fa-a291-d66ebca29c21/scratchpad/qllm")
from model import load, Tokenizer, forward_fp

BITS_PER_U = 16
PROMPT = "Once upon a time, there was a little"
TEMP = 1.0
TOP_P = 0.9
SEED_MT = 20261108
SEED_CS = 7


class Pool:
    """按 16 bit 产出 [0,1) 均匀数的随机源。"""

    def __init__(self, kind, bits=None):
        self.kind = kind
        self.i = 0
        self.n_used = 0
        if kind == "QRNG":
            self.bits = bits
            self.cap = bits.size // BITS_PER_U
        elif kind == "MT19937":
            self.rng = random.Random(SEED_MT)
        elif kind == "CSPRNG":
            pass

    def u(self):
        self.n_used += 1
        if self.kind == "QRNG":
            if self.i + BITS_PER_U > self.bits.size:
                raise RuntimeError("QRNG pool exhausted")
            chunk = self.bits[self.i:self.i + BITS_PER_U]
            self.i += BITS_PER_U
            v = int("".join(map(str, chunk)), 2)
        elif self.kind == "MT19937":
            v = self.rng.getrandbits(BITS_PER_U)
        else:  # CSPRNG
            v = int.from_bytes(os.urandom(2), "big")
        return (v + 0.5) / (1 << BITS_PER_U)

    def remaining(self):
        if self.kind == "QRNG":
            return self.cap - self.i // BITS_PER_U
        return float("inf")


def topp_dist(logits, temp=TEMP, top_p=TOP_P):
    z = logits / temp
    z = z - z.max()
    p = np.exp(z)
    p /= p.sum()
    order = np.argsort(p)[::-1]
    cum = np.cumsum(p[order])
    k = int(np.searchsorted(cum, top_p)) + 1
    keep = order[:k]
    q = np.zeros_like(p)
    q[keep] = p[keep]
    q /= q.sum()
    return q, keep


def sample_from(dist, u):
    return int(np.searchsorted(np.cumsum(dist), u, side="right"))


def experiment_A(cfg, W, tok, pools, N):
    ids = tok.encode(PROMPT)
    logits = forward_fp(cfg, W, ids)[-1]
    P, keep = topp_dist(logits)
    support = np.sort(keep)
    Ps = P[support]
    res = {}
    for name, pool in pools.items():
        counts = np.zeros(support.size, dtype=float)
        idx_of = {int(t): i for i, t in enumerate(support)}
        us = []
        for _ in range(N):
            u = pool.u()
            us.append(u)
            t = sample_from(P, u)
            counts[idx_of[t]] += 1
        exp = N * Ps
        # 合并期望计数过小的类别，保证卡方有效
        order = np.argsort(exp)[::-1]
        oc, oe = counts[order], exp[order]
        mc, me = [], []
        acc_c = acc_e = 0.0
        for c, e in zip(oc, oe):
            acc_c += c
            acc_e += e
            if acc_e >= 5:
                mc.append(acc_c)
                me.append(acc_e)
                acc_c = acc_e = 0.0
        if acc_e > 0:
            mc[-1] += acc_c
            me[-1] += acc_e
        me = np.array(me) * (sum(mc) / sum(me))
        chi2_p = chisquare(mc, me).pvalue
        ks_p = kstest(us, "uniform").pvalue
        tv = 0.5 * np.sum(np.abs(counts / N - Ps))
        res[name] = dict(support=int(support.size), draws=N, chi2_gof_p=float(chi2_p),
                         ks_uniform_p=float(ks_p), tv_distance=float(tv),
                         u_mean=float(np.mean(us)), u_std=float(np.std(us)))
    return res, dict(support_tokens=[tok.decode([int(t)]) for t in support],
                     target_probs=[float(x) for x in Ps])


def generate(cfg, W, tok, pool, n_cont, n_new):
    base = tok.encode(PROMPT)
    conts = []
    for _ in range(n_cont):
        seq = list(base)
        for _ in range(n_new):
            logits = forward_fp(cfg, W, seq)[-1]
            P, _ = topp_dist(logits)
            seq.append(sample_from(P, pool.u()))
        conts.append(seq[len(base):])
    return conts


def distinct_n(conts, n):
    grams = {}
    tot = 0
    for c in conts:
        for i in range(len(c) - n + 1):
            grams[tuple(c[i:i + n])] = 1
            tot += 1
    return len(grams) / max(tot, 1)


def experiment_B(cfg, W, tok, pools, n_cont, n_new):
    res = {}
    all_tokens = {}
    for name, pool in pools.items():
        conts = generate(cfg, W, tok, pool, n_cont, n_new)
        flat = [t for c in conts for t in c]
        uniq = len(set(tuple(c) for c in conts)) / len(conts)
        tf = np.bincount(flat, minlength=cfg.vocab).astype(float)
        all_tokens[name] = tf
        res[name] = dict(n_cont=n_cont, n_new=n_new,
                         distinct_1=distinct_n(conts, 1), distinct_2=distinct_n(conts, 2),
                         distinct_3=distinct_n(conts, 3), unique_frac=uniq,
                         tokens_generated=len(flat), u_used=pool.n_used)
        res[name]["sample_texts"] = [tok.decode(c) for c in conts[:3]]
    # 源间 token 频率相关性
    names = list(all_tokens)
    corr = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = all_tokens[names[i]], all_tokens[names[j]]
            m = (a + b) > 0
            corr[f"{names[i]}~{names[j]}"] = float(np.corrcoef(a[m], b[m])[0, 1])
    return res, corr


def main():
    cfg, W = load()
    tok = Tokenizer()
    ext = np.load("extracted_bits.npy")
    print(f"QRNG pool: {ext.size} bit -> {ext.size // BITS_PER_U} uniforms (16-bit)")

    # 预算：实验 A 15000 draws + 实验 B 500×20=10000，QRNG 共需 25000 < 34350
    NA = 15000
    NB_cont, NB_new = 500, 20

    poolsA = {"QRNG": Pool("QRNG", ext), "MT19937": Pool("MT19937"), "CSPRNG": Pool("CSPRNG")}
    resA, meta = experiment_A(cfg, W, tok, poolsA, NA)
    print("\n实验 A：首 token 分布保真度（卡方/KS p>=0.01 为通过）")
    for k, v in resA.items():
        print(f"  {k:9s} chi2_gof_p={v['chi2_gof_p']:.3f} ks_uniform_p={v['ks_uniform_p']:.3f} TV={v['tv_distance']:.4f}")

    # 实验 B 用新的 QRNG 池游标（接着 A 继续消耗同一真机池）
    poolsB = {"QRNG": poolsA["QRNG"], "MT19937": Pool("MT19937"), "CSPRNG": Pool("CSPRNG")}
    print(f"\nQRNG remaining before B: {poolsB['QRNG'].remaining()} uniforms")
    resB, corr = experiment_B(cfg, W, tok, poolsB, NB_cont, NB_new)
    print("\n实验 B：多样性")
    for k, v in resB.items():
        print(f"  {k:9s} distinct1={v['distinct_1']:.3f} distinct2={v['distinct_2']:.3f} "
              f"distinct3={v['distinct_3']:.3f} unique={v['unique_frac']:.3f}")
    print("  token频率相关:", {k: round(v, 4) for k, v in corr.items()})

    out = dict(prompt=PROMPT, temp=TEMP, top_p=TOP_P, bits_per_uniform=BITS_PER_U,
               qrng_pool_uniforms=int(ext.size // BITS_PER_U),
               experiment_A=resA, experiment_A_meta=meta,
               experiment_B=resB, token_freq_corr=corr)
    json.dump(out, open("qsample_results.json", "w"), ensure_ascii=False, indent=1)
    print("\nwrote qsample_results.json")


if __name__ == "__main__":
    main()
