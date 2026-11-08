"""NIST SP 800-22 随机性检验子集（从零实现，仅依赖 numpy / scipy.special）。

每个检验返回 p 值；按惯例 p >= 0.01 视为"未检出非随机性"（通过）。
实现覆盖：单比特频数、块内频数、游程、块内最长游程、二元矩阵秩、
离散傅里叶（谱）、累加和、近似熵、序列（Serial）。
参考：NIST SP 800-22 Rev.1a。
"""
import numpy as np
from scipy.special import erfc, gammaincc, gammainc
from math import sqrt, log, floor


def _bits(x):
    return np.asarray(x, dtype=np.int8)


def frequency_monobit(bits):
    b = _bits(bits)
    n = b.size
    s = np.sum(2 * b - 1)
    sobs = abs(s) / sqrt(n)
    return erfc(sobs / sqrt(2))


def block_frequency(bits, M=128):
    b = _bits(bits)
    n = b.size
    N = n // M
    if N == 0:
        return float('nan')
    blocks = b[:N * M].reshape(N, M)
    pi = blocks.mean(axis=1)
    chi2 = 4.0 * M * np.sum((pi - 0.5) ** 2)
    return gammaincc(N / 2.0, chi2 / 2.0)


def runs(bits):
    b = _bits(bits)
    n = b.size
    pi = b.mean()
    if abs(pi - 0.5) >= 2.0 / sqrt(n):
        return 0.0  # 先验频数检验未过，游程检验不适用
    vobs = 1 + np.sum(b[1:] != b[:-1])
    num = abs(vobs - 2 * n * pi * (1 - pi))
    den = 2 * sqrt(2 * n) * pi * (1 - pi)
    return erfc(num / den)


def longest_run_of_ones(bits):
    b = _bits(bits)
    n = b.size
    # 按 n 选块参数（SP800-22 表）
    if n < 128:
        return float('nan')
    if n < 6272:
        M, K, Nblk = 8, 3, 16
        v_lo, pi = 1, [0.2148, 0.3672, 0.2305, 0.1875]
    elif n < 750000:
        M, K, Nblk = 128, 5, 49
        v_lo, pi = 4, [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
    else:
        M, K, Nblk = 10000, 6, 75
        v_lo, pi = 10, [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    blocks = b[:Nblk * M].reshape(Nblk, M)
    counts = np.zeros(K + 1, dtype=int)
    for blk in blocks:
        mx = 0
        cur = 0
        for bit in blk:
            if bit == 1:
                cur += 1
                mx = max(mx, cur)
            else:
                cur = 0
        idx = min(max(mx, v_lo), v_lo + K) - v_lo
        counts[idx] += 1
    pi = np.array(pi)
    chi2 = np.sum((counts - Nblk * pi) ** 2 / (Nblk * pi))
    return gammaincc(K / 2.0, chi2 / 2.0)


def binary_matrix_rank(bits, M=32, Q=32):
    b = _bits(bits)
    n = b.size
    N = n // (M * Q)
    if N == 0:
        return float('nan')

    def rank_gf2(mat):
        m = mat.copy().astype(np.int8)
        rows, cols = m.shape
        r = 0
        for c in range(cols):
            piv = None
            for i in range(r, rows):
                if m[i, c]:
                    piv = i
                    break
            if piv is None:
                continue
            m[[r, piv]] = m[[piv, r]]
            for i in range(rows):
                if i != r and m[i, c]:
                    m[i] ^= m[r]
            r += 1
            if r == rows:
                break
        return r

    FM = FM1 = 0
    rem = 0
    for k in range(N):
        mat = b[k * M * Q:(k + 1) * M * Q].reshape(M, Q)
        rk = rank_gf2(mat)
        if rk == M:
            FM += 1
        elif rk == M - 1:
            FM1 += 1
        else:
            rem += 1
    p_full, p_f1 = 0.2888, 0.5776
    p_rem = 1 - p_full - p_f1
    chi2 = ((FM - p_full * N) ** 2 / (p_full * N)
            + (FM1 - p_f1 * N) ** 2 / (p_f1 * N)
            + (rem - p_rem * N) ** 2 / (p_rem * N))
    return np.exp(-chi2 / 2.0)


def dft_spectral(bits):
    b = _bits(bits)
    n = b.size
    x = 2 * b - 1
    s = np.abs(np.fft.fft(x))[:n // 2]
    T = sqrt(log(1.0 / 0.05) * n)
    N0 = 0.95 * n / 2.0
    N1 = np.sum(s < T)
    d = (N1 - N0) / sqrt(n * 0.95 * 0.05 / 4.0)
    return erfc(abs(d) / sqrt(2))


def cumulative_sums(bits, mode=0):
    b = _bits(bits)
    n = b.size
    x = 2 * b - 1
    if mode == 0:
        z = np.max(np.abs(np.cumsum(x)))
    else:
        z = np.max(np.abs(np.cumsum(x[::-1])))
    if z == 0:
        return 1.0
    from scipy.stats import norm
    start = int(floor((-n / z + 1) / 4.0))
    end = int(floor((n / z - 1) / 4.0))
    s1 = sum(norm.cdf((4 * k + 1) * z / sqrt(n)) - norm.cdf((4 * k - 1) * z / sqrt(n))
             for k in range(start, end + 1))
    start = int(floor((-n / z - 3) / 4.0))
    s2 = sum(norm.cdf((4 * k + 3) * z / sqrt(n)) - norm.cdf((4 * k + 1) * z / sqrt(n))
             for k in range(start, end + 1))
    return 1.0 - s1 + s2


def _psi2(b, m):
    n = b.size
    if m == 0:
        return 0.0
    ext = np.concatenate([b, b[:m - 1]])
    # 把每个长 m 窗口当作二进制数计数
    idx = np.zeros(n, dtype=np.int64)
    for j in range(m):
        idx = (idx << 1) | ext[j:j + n]
    counts = np.bincount(idx, minlength=1 << m).astype(float)
    return np.sum(counts ** 2) * (1 << m) / n - n


def serial(bits, m=16):
    b = _bits(bits)
    p0 = _psi2(b, m)
    p1 = _psi2(b, m - 1)
    p2 = _psi2(b, m - 2)
    d1 = p0 - p1
    d2 = p0 - 2 * p1 + p2
    pv1 = gammaincc(2 ** (m - 2), d1 / 2.0)
    pv2 = gammaincc(2 ** (m - 3), d2 / 2.0)
    return pv1, pv2


def approximate_entropy(bits, m=10):
    b = _bits(bits)
    n = b.size

    def phi(mm):
        if mm == 0:
            return 0.0
        ext = np.concatenate([b, b[:mm - 1]])
        idx = np.zeros(n, dtype=np.int64)
        for j in range(mm):
            idx = (idx << 1) | ext[j:j + n]
        counts = np.bincount(idx, minlength=1 << mm).astype(float)
        c = counts / n
        nz = c[c > 0]
        return np.sum(nz * np.log(nz))

    apen = phi(m) - phi(m + 1)
    chi2 = 2.0 * n * (log(2) - apen)
    return gammaincc(2 ** (m - 1), chi2 / 2.0)


def run_all(bits, label=""):
    b = _bits(bits)
    out = {}
    out["frequency_monobit"] = frequency_monobit(b)
    out["block_frequency_M128"] = block_frequency(b, 128)
    out["runs"] = runs(b)
    out["longest_run_ones"] = longest_run_of_ones(b)
    out["binary_matrix_rank"] = binary_matrix_rank(b)
    out["dft_spectral"] = dft_spectral(b)
    out["cusum_forward"] = cumulative_sums(b, 0)
    out["cusum_backward"] = cumulative_sums(b, 1)
    s1, s2 = serial(b, 16)
    out["serial_1"] = s1
    out["serial_2"] = s2
    out["approx_entropy_m10"] = approximate_entropy(b, 10)
    return out
