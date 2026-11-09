"""任务 #22：对 QZ03 真机原始比特做偏置/最小熵估计、Toeplitz 去偏提取、NIST 检验（原始 vs 提取）。

- 原始流：20000 shots × 35 qubit，按 shot 行优先展平 = 700000 bit。
- 最小熵：逐比特 h_min = -log2(max(p0,p1))，用 99% 置信上界对 p_max 做保守估计。
- 提取：按块 Toeplitz 哈希（强提取器），m×n 矩阵由 m+n-1 个种子比特定义（种子取自 CSPRNG，公开）。
  提取率 m/n 取保守最小熵率减安全余量（leftover hash lemma，ε=2^-100 → 损失 200 bit）。
- 输出：qrng_analysis.json（全部数值），extracted_bits.npy（提取后比特，供 #23/#24 使用），raw_bits.npy。
"""
import json
import numpy as np
from math import log2, sqrt
from scipy.stats import norm

import nist_sts as nist

RAW = "qz03_qrng_raw.json"


def load_raw():
    d = json.load(open(RAW))
    words = [s for t in d["tasks"] for s in t["memory"]]
    B = np.array([[int(c) for c in s] for s in words], dtype=np.int8)  # [shots, 35]
    return d, B


def min_entropy(B, conf=0.99):
    """逐比特最小熵：经验值 + 保守（置信上界 p_max）值。"""
    N = B.shape[0]
    z = norm.ppf(conf)
    p1 = B.mean(axis=0)
    pmax = np.maximum(p1, 1 - p1)
    # 二项比例置信上界，保守抬高 pmax -> 压低最小熵
    pmax_u = np.minimum(1.0, pmax + z * np.sqrt(pmax * (1 - pmax) / N))
    h_emp = -np.log2(pmax)
    h_cons = -np.log2(pmax_u)
    return dict(p1=p1, pmax=pmax, pmax_u=pmax_u, h_emp=h_emp, h_cons=h_cons,
                rate_emp=float(h_emp.mean()), rate_cons=float(h_cons.mean()),
                total_emp=float(h_emp.sum()), total_cons=float(h_cons.sum()))


def toeplitz_extract(stream, n_in, m_out, seed=1234567):
    """按块 Toeplitz 提取。矩阵 m_out×n_in，由 m_out+n_in-1 个种子比特定义；块间共用同一矩阵。"""
    rng = np.random.default_rng(seed)
    tcol = rng.integers(0, 2, size=m_out + n_in - 1, dtype=np.int8)  # Toeplitz 由一条对角线定义
    nblocks = stream.size // n_in
    S = stream[:nblocks * n_in].reshape(nblocks, n_in).astype(np.int8)
    # 构造 Toeplitz 矩阵 T[i,j] = tcol[i - j + (n_in-1)]
    i = np.arange(m_out)[:, None]
    j = np.arange(n_in)[None, :]
    T = tcol[i - j + (n_in - 1)]  # [m_out, n_in]
    # 每块输出 = T · block (mod 2)
    out = (S @ T.T) & 1  # [nblocks, m_out]
    return out.reshape(-1).astype(np.int8), dict(n_in=n_in, m_out=m_out, nblocks=nblocks,
                                                 seed=seed, out_bits=int(out.size))


def main():
    d, B = load_raw()
    shots, nq = B.shape
    raw = B.reshape(-1).astype(np.int8)  # 行优先：shot0 的 35 位, shot1 的 35 位, ...
    np.save("raw_bits.npy", raw)

    me = min_entropy(B)
    print(f"shots={shots} qubits={nq} raw_bits={raw.size}")
    print(f"ones_frac={raw.mean():.4f}")
    print(f"min-entropy rate: empirical={me['rate_emp']:.4f} conservative(99%)={me['rate_cons']:.4f} bit/bit")
    print(f"total conservative min-entropy ~ {me['total_cons']*shots/1:.0f}? per word {me['total_cons']:.2f} bit")

    # 提取参数：输入块 = 一个整段 2800 bit（= 80 shots × 35），输出率 < 保守最小熵率 - 余量
    eps_loss_bits = 200  # 2 * log2(1/eps), eps = 2^-100
    n_in = 35 * 200          # 7000 bit 一块
    budget = me["rate_cons"] * n_in - eps_loss_bits
    m_out = int(budget * 0.95)  # 再留 5% 工程余量
    m_out -= m_out % 8
    ext, einfo = toeplitz_extract(raw, n_in, m_out, seed=20261108)
    np.save("extracted_bits.npy", ext)
    print(f"extract: block {n_in}->{m_out} bit (rate {m_out/n_in:.3f}), blocks={einfo['nblocks']}, out={ext.size} bit")
    print(f"extracted ones_frac={ext.mean():.5f}")

    # NIST 检验：原始 vs 提取（取等长前缀便于对照）
    L = min(raw.size, ext.size)
    L -= L % (32 * 32)   # 便于矩阵秩检验分块
    nist_raw = nist.run_all(raw[:L])
    nist_ext = nist.run_all(ext[:L])
    print("\nNIST SP 800-22 (p>=0.01 通过)")
    print(f"{'test':24s} {'raw':>10s} {'extracted':>12s}")
    for k in nist_raw:
        print(f"{k:24s} {nist_raw[k]:10.4f} {nist_ext[k]:12.4f}")

    raw_pass = sum(1 for v in nist_raw.values() if v >= 0.01)
    ext_pass = sum(1 for v in nist_ext.values() if v >= 0.01)
    print(f"\nraw passes {raw_pass}/{len(nist_raw)}, extracted passes {ext_pass}/{len(nist_ext)}")

    out = dict(
        shots=shots, qubits=nq, raw_bits=int(raw.size),
        total_qpu_exec_s=float(sum(t["execution_time"] for t in d["tasks"])),
        total_cost=float(sum(t["cost"] for t in d["tasks"])),
        ones_frac_raw=float(raw.mean()), ones_frac_ext=float(ext.mean()),
        per_qubit_p1=me["p1"].round(5).tolist(),
        per_qubit_hmin_emp=me["h_emp"].round(5).tolist(),
        per_qubit_hmin_cons=me["h_cons"].round(5).tolist(),
        min_entropy_rate_emp=me["rate_emp"], min_entropy_rate_cons=me["rate_cons"],
        min_entropy_per_word_cons=me["total_cons"],
        extraction=einfo, extraction_rate=float(m_out / n_in),
        nist_raw=nist_raw, nist_ext=nist_ext,
        nist_raw_pass=raw_pass, nist_ext_pass=ext_pass, nist_total=len(nist_raw),
        tests_len=int(L),
    )
    json.dump(out, open("qrng_analysis.json", "w"), ensure_ascii=False, indent=1)
    print("\nwrote qrng_analysis.json, raw_bits.npy, extracted_bits.npy")


if __name__ == "__main__":
    main()
