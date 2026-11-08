"""任务 #25：吞吐、延迟与成本分析。

数据来源：qz03_qrng_raw.json（20 条真机任务的 created/started/completed 时间戳、execution_time、cost）。
量子侧每条 = 35 qubit × 1000 shots = 35000 原始 bit；经 Toeplitz 提取率 0.785 → ~27475 可用 bit。
分析：
  - 原始/提取比特率（按 QPU 纯执行时间，和按含排队的墙钟时间两种口径）
  - 单比特成本
  - 每个采样 token 的随机数开销与等效延迟（池化后，采样本身只是查表，QPU 仅负责离线补充池）
  - 与光子 QRNG / PRNG 的数量级对照（引用正文参考文献，不臆造精确值）
"""
import json
import numpy as np
from datetime import datetime

RAW = "qz03_qrng_raw.json"
EXTRACT_RATE = 0.785
BITS_PER_TASK = 35 * 1000
BITS_PER_U = 16


def parse(ts):
    return datetime.fromisoformat(ts.replace("Z", ""))


def main():
    d = json.load(open(RAW))
    T = d["tasks"]
    exec_s = np.array([t["execution_time"] for t in T])
    cost = np.array([t["cost"] for t in T])

    # 墙钟延迟（created -> completed），含排队与网络
    wall = []
    for t in T:
        try:
            wall.append((parse(t["completed_at"]) - parse(t["created_at"])).total_seconds())
        except Exception:
            pass
    wall = np.array(wall, dtype=float)

    n = len(T)
    raw_bits = n * BITS_PER_TASK
    ext_bits = int(raw_bits * EXTRACT_RATE)

    # 吞吐：纯 QPU 执行口径
    rate_exec_raw = BITS_PER_TASK / exec_s.mean()          # bit/s（单条）
    rate_exec_ext = rate_exec_raw * EXTRACT_RATE
    # 吞吐：含排队墙钟口径
    rate_wall_raw = BITS_PER_TASK / wall.mean() if wall.size else float("nan")

    cost_per_raw_bit = cost.sum() / raw_bits               # 元/bit
    cost_per_ext_bit = cost.sum() / ext_bits
    cost_per_megabit_ext = cost_per_ext_bit * 1e6

    # 采样开销：每个 token 采样消耗 1 个 16-bit 均匀数
    bits_per_token = BITS_PER_U
    ext_bits_total = ext_bits
    tokens_supported = ext_bits_total // bits_per_token
    # 池化后采样延迟 = 查表（纳秒级），QPU 只在后台离线补池；
    # 维持某 token 生成速率所需的 QPU 比特产率：
    def qpu_rate_needed(tok_per_s):
        return tok_per_s * bits_per_token / EXTRACT_RATE   # 需要的原始 bit/s

    out = dict(
        n_tasks=n,
        raw_bits_total=raw_bits, extracted_bits_total=ext_bits,
        exec_s_mean=float(exec_s.mean()), exec_s_std=float(exec_s.std()),
        exec_s_total=float(exec_s.sum()),
        wall_s_mean=float(wall.mean()) if wall.size else None,
        wall_s_min=float(wall.min()) if wall.size else None,
        wall_s_max=float(wall.max()) if wall.size else None,
        rate_exec_raw_bps=float(rate_exec_raw), rate_exec_ext_bps=float(rate_exec_ext),
        rate_wall_raw_bps=float(rate_wall_raw),
        cost_total=float(cost.sum()), currency=T[0]["currency"],
        cost_per_raw_bit=float(cost_per_raw_bit), cost_per_ext_bit=float(cost_per_ext_bit),
        cost_per_megabit_ext=float(cost_per_megabit_ext),
        bits_per_token=bits_per_token, tokens_supported_by_pool=int(tokens_supported),
        qpu_raw_bps_needed_for_1ktok_s=float(qpu_rate_needed(1000)),
        qpu_raw_bps_needed_for_1tok_s=float(qpu_rate_needed(1)),
    )
    print(json.dumps(out, ensure_ascii=False, indent=1))
    json.dump(out, open("throughput_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n--- 关键数字 ---")
    print(f"每条任务 QPU 执行 {out['exec_s_mean']:.2f}±{out['exec_s_std']:.2f} s，产 {BITS_PER_TASK} 原始 bit")
    print(f"墙钟延迟（含排队）{out['wall_s_mean']:.1f} s（{out['wall_s_min']:.0f}–{out['wall_s_max']:.0f}）")
    print(f"纯执行吞吐：原始 {rate_exec_raw:,.0f} bit/s，提取后 {rate_exec_ext:,.0f} bit/s")
    print(f"成本：¥{out['cost_total']:.0f} 共 {ext_bits:,} 提取 bit → ¥{cost_per_megabit_ext:.2f}/Mbit")
    print(f"池可支撑 {tokens_supported:,} 次 token 采样（16 bit/token）")
    print(f"维持 1000 token/s 需 QPU 原始产率 {qpu_rate_needed(1000):,.0f} bit/s")


if __name__ == "__main__":
    main()
