"""论文配图（matplotlib，Noto Sans CJK SC，色盲友好配色）。"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
fm.fontManager.addfont(FONT)
plt.rcParams["font.family"] = fm.FontProperties(fname=FONT).get_name()
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

# 色盲友好
C = dict(blue="#4E79A7", orange="#F28E2B", green="#59A14F", red="#E15759",
         gray="#9C9C9C", purple="#B07AA1")

A = json.load(open("qrng_analysis.json"))
TH = json.load(open("throughput_results.json"))


def fig_bias():
    p1 = np.array(A["per_qubit_p1"])
    hcons = np.array(A["per_qubit_hmin_cons"])
    q = np.arange(len(p1))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
    ax[0].bar(q, p1, color=[C["red"] if abs(x - 0.5) > 0.05 else C["blue"] for x in p1])
    ax[0].axhline(0.5, color=C["gray"], lw=1, ls="--")
    ax[0].set_xlabel("逻辑比特编号")
    ax[0].set_ylabel("输出为 1 的频率 $p_1$")
    ax[0].set_title("(a) QZ03 各比特偏置（20000 shots）")
    ax[0].set_ylim(0, 0.6)
    ax[1].bar(q, hcons, color=C["blue"])
    ax[1].axhline(hcons.mean(), color=C["orange"], lw=1.5,
                  label=f"平均 {hcons.mean():.3f} bit")
    ax[1].set_xlabel("逻辑比特编号")
    ax[1].set_ylabel("最小熵 $H_\\infty$ (bit)")
    ax[1].set_title("(b) 逐比特最小熵（99% 保守下界）")
    ax[1].set_ylim(0, 1.05)
    ax[1].legend(frameon=False, fontsize=9)
    fig.savefig("fig_bias.pdf")
    plt.close(fig)


def fig_nist():
    raw = A["nist_raw"]
    ext = A["nist_ext"]
    names = list(raw.keys())
    labels = {"frequency_monobit": "单比特频数", "block_frequency_M128": "块内频数",
              "runs": "游程", "longest_run_ones": "最长游程", "binary_matrix_rank": "矩阵秩",
              "dft_spectral": "谱(DFT)", "cusum_forward": "累加和+", "cusum_backward": "累加和-",
              "serial_1": "序列1", "serial_2": "序列2", "approx_entropy_m10": "近似熵"}
    y = np.arange(len(names))
    r = [max(raw[n], 1e-4) for n in names]
    e = [max(ext[n], 1e-4) for n in names]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    h = 0.38
    ax.barh(y + h / 2, r, h, color=C["red"], label="原始真机比特")
    ax.barh(y - h / 2, e, h, color=C["green"], label="Toeplitz 提取后")
    ax.axvline(0.01, color=C["gray"], ls="--", lw=1.2)
    ax.text(0.0115, len(names) - 0.5, "通过阈值 0.01", fontsize=8, color=C["gray"])
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([labels[n] for n in names], fontsize=9)
    ax.set_xlabel("NIST SP 800-22  $p$ 值（对数轴，越大越随机）")
    ax.set_title(f"随机性检验：原始 {A['nist_raw_pass']}/{A['nist_total']} 通过 → "
                 f"提取后 {A['nist_ext_pass']}/{A['nist_total']} 通过")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    ax.set_xlim(8e-5, 1.5)
    fig.savefig("fig_nist.pdf")
    plt.close(fig)


def fig_fidelity():
    S = json.load(open("qsample_results.json")) if _exists("qsample_results.json") else None
    if S is None:
        return
    meta = S["experiment_A_meta"]
    toks = meta["support_tokens"]
    P = np.array(meta["target_probs"])
    order = np.argsort(P)[::-1][:12]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    x = np.arange(len(order))
    ax.bar(x, P[order], color=C["gray"], width=0.6, label="目标分布 $P$（top-p）", zorder=1)
    resA = S["experiment_A"]
    # 叠加三源经验频率点（用 TV 已知一致，这里画 chi2/ks 的说明在表里；点示意）
    ax.set_xticks(x)
    ax.set_xticklabels([repr(toks[i])[1:-1][:8] for i in order], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("概率")
    ax.set_title("首 token 的 top-p 目标分布（三种随机源均复现，见表）")
    ax.legend(frameon=False, fontsize=9)
    fig.savefig("fig_fidelity.pdf")
    plt.close(fig)


def fig_attack():
    AT = json.load(open("attack_results.json"))
    mt_num = AT["mt_number_prediction"]["hit_rate"] * 100
    mt_tok = AT["mt_token_prediction"]["hit_rate"] * 100
    qr_tok = AT["qrng_token_prediction"]["hit_rate"] * 100
    base = AT["qrng_token_prediction"]["collision_baseline"] * 100
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    cats = ["MT19937\n数值预测", "MT19937\ntoken 预测", "QRNG 池\ntoken 预测"]
    vals = [mt_num, mt_tok, qr_tok]
    cols = [C["red"], C["red"], C["green"]]
    bars = ax.bar(cats, vals, color=cols, width=0.6)
    ax.axhline(base, color=C["gray"], ls="--", lw=1.3)
    ax.text(2.05, base + 1, f"盲猜碰撞基线 {base:.0f}%", fontsize=8, color=C["gray"], ha="right")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}%", ha="center", fontsize=10)
    ax.set_ylabel("攻击者命中率")
    ax.set_ylim(0, 110)
    ax.set_title("种子劫持/预测攻击：观测 624 个输出后")
    fig.savefig("fig_attack.pdf")
    plt.close(fig)


def _exists(p):
    import os
    return os.path.exists(p)


if __name__ == "__main__":
    fig_bias()
    fig_nist()
    fig_attack()
    fig_fidelity()
    print("figures written (bias, nist, attack" + (", fidelity" if _exists("qsample_results.json") else "") + ")")
