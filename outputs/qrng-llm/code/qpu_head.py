"""在逻辑比特云 QZ03 真机上计算 stories15M 第 0 层第 0 个注意力头的打分。

场景：模型写到 "Once upon a time, there was a little ___"，"little" 这个位置的查询 q
对前面 10 个位置的键 k_j 打分。每个 48 维内积拆成 12 段 4 维内积：
    q·k_j = Σ_c |q_c| |k_jc| <û_c | v̂_jc>
每段 <û|v̂> 用一条 3 比特干涉线路测：(|0>|û> + |1>|v̂>)/√2 → 辅助比特 H → 测辅助比特，
P0 - P1 = <û|v̂>。一条大线路在 10 组互不相交的相邻比特上并行跑 10 段。
校准线路：已知 <û|v̂> = +1 / -1 的段，用来拟合每组比特的对比度与偏置。

本文件只在本地构造线路，并用 SDK 的官方序列化器生成提交用的 JSON（不联网、不含 API key）。
"""
import json
import math
import sys
import numpy as np
from qiskit import QuantumCircuit as QQC, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector

import lqcloud.provider.cloud_provider as cp
from lqcloud import LQCloudProvider, QuantumCircuit as LQC

from model import load, Tokenizer, rmsnorm, rope_tables, apply_rope

QZ03 = {"name": "QZ03", "type": "real_qpu", "status": "active", "qubits": 35, "native_gates": None,
        "supported_modes": ["circuit", "wave", "sequence"], "supported_measurement_types": ["measure2"],
        "topology": {"type": "grid", "rows": 6, "cols": 6, "coupling_map": [
            [0, 1], [0, 6], [1, 2], [1, 7], [2, 3], [2, 8], [3, 4], [3, 9], [4, 10], [5, 11], [6, 7], [6, 12], [7, 8], [7, 13],
            [8, 9], [9, 10], [10, 11], [10, 15], [11, 16], [12, 13], [12, 17], [14, 15], [14, 20], [15, 16], [16, 22], [17, 18],
            [17, 23], [18, 19], [18, 24], [19, 20], [20, 26], [21, 22], [21, 27], [22, 28], [23, 24], [23, 29], [24, 25],
            [24, 30], [25, 26], [25, 31], [26, 27], [26, 32], [27, 33], [28, 34], [29, 30], [30, 31], [31, 32], [32, 33]]}}
# 10 组相邻比特链 [辅助, 数据1, 数据2]，避开保真度低的耦合
TRIPLES = [[0, 1, 2], [6, 7, 13], [3, 4, 10], [5, 11, 16], [12, 17, 18],
           [19, 20, 14], [23, 24, 30], [32, 31, 25], [26, 27, 33], [21, 22, 28]]
BLOCK = 4


def block_circuit(u, v):
    """3 比特 Hadamard 测试，返回转译到 {rz, sx, x, cz}、线性连接的 Qiskit 线路（不含测量）。"""
    psi = np.zeros(8)
    for d in range(4):
        psi[2 * d] = u[d] / math.sqrt(2)
        psi[2 * d + 1] = v[d] / math.sqrt(2)
    qc = QQC(3)
    qc.append(StatePreparation(psi), [0, 1, 2])
    qc.h(0)
    t = transpile(qc, basis_gates=["rz", "sx", "x", "cz"], coupling_map=[[0, 1], [1, 0], [1, 2], [2, 1]],
                  initial_layout=[0, 1, 2], layout_method="trivial", optimization_level=3, seed_transpiler=7)
    fl = t.layout.final_index_layout() if t.layout is not None else [0, 1, 2]
    anc = fl[0]                      # 路由后辅助比特最终所在的位置
    p = Statevector(t).probabilities([anc])
    assert abs((p[0] - p[1]) - float(u @ v)) < 1e-9
    return t, anc


def model_case():
    cfg, W = load()
    tok = Tokenizer()
    ids = tok.encode("Once upon a time, there was a little")
    T = len(ids)
    xn = rmsnorm(W["tok"][ids], W["rms_att"][0])
    cos, sin = rope_tables(cfg, T)
    q = apply_rope(xn @ W["wq"][0].T, cos, sin)[-1, :48]
    K = apply_rope(xn @ W["wk"][0].T, cos, sin)[:, :48]
    pieces = [tok.vocab[t].decode("utf-8", "replace") if t != 1 else "<s>" for t in ids]
    return q, K, pieces


def build_tests(q, K):
    tests = []   # (key j, block c, weight, û, v̂)
    for j in range(len(K)):
        for c in range(48 // BLOCK):
            qc_, kc = q[c * BLOCK:(c + 1) * BLOCK], K[j, c * BLOCK:(c + 1) * BLOCK]
            tests.append(dict(j=j, c=c, w=float(np.linalg.norm(qc_) * np.linalg.norm(kc)),
                              u=qc_ / np.linalg.norm(qc_), v=kc / np.linalg.norm(kc)))
    return tests


def pack(tests_10):
    """把至多 10 段测试打包成一条 30 逻辑比特的 lqcloud 线路，辅助比特 t 测到经典位 t。"""
    lq = LQC(3 * len(TRIPLES), len(TRIPLES))
    n_cz = 0
    ancs = []
    for t_i, tst in enumerate(tests_10):
        tq, anc = block_circuit(tst["u"], tst["v"])
        base = 3 * t_i
        ancs.append(base + anc)
        for inst in tq.data:
            name = inst.operation.name
            qs = [base + tq.find_bit(b).index for b in inst.qubits]
            if name == "rz":
                lq.rz(float(inst.operation.params[0]), qs[0])
            elif name == "sx":
                lq.sx(qs[0])
            elif name == "x":
                lq.x(qs[0])
            elif name == "cz":
                lq.cz(qs[0], qs[1]); n_cz += 1
            else:
                raise ValueError(name)
    lq.barrier()
    for t_i, a in enumerate(ancs):
        lq.measure(a, t_i)
    return lq, n_cz


def capture_payload(circuits, shots):
    """用 SDK 的官方序列化流程生成提交 JSON，拦截 HTTP，不真正发送。"""
    box = {}

    class _Stop(Exception):
        pass

    def fake(session, method, url, **kw):
        box.update(method=method, url=url, payload=kw.get("json"))
        raise _Stop()

    cp.request_with_retry = fake
    prov = LQCloudProvider(api_key="offline-placeholder", url="https://cloud.logicalqubit.com", interactive=False)
    prov.get_backend_config = lambda name: QZ03
    be = prov.get_backend("QZ03", verify=False)
    be._config_cache = QZ03
    layout = [p for tr in TRIPLES for p in tr]
    try:
        be.run(circuits if len(circuits) > 1 else circuits[0], shots=shots, initial_layout=layout)
    except _Stop:
        pass
    return box


def calibration_tests(u_ref, sign_even):
    out = []
    for t_i in range(len(TRIPLES)):
        s = sign_even if t_i % 2 == 0 else -sign_even
        out.append(dict(j=-1, c=-1, w=0.0, u=u_ref, v=s * u_ref, expect=float(s)))
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    q, K, pieces = model_case()
    tests = build_tests(q, K)
    u_ref = tests[0]["u"]
    cal_a, cal_b = calibration_tests(u_ref, +1.0), calibration_tests(u_ref, -1.0)
    if mode == "smoke":
        lq, ncz = pack(cal_a)
        box = capture_payload([lq], shots=500)
        json.dump(dict(plan=[dict(kind="cal", expect=t["expect"]) for t in cal_a], **box),
                  open("qpu_smoke_payload.json", "w"))
        print("smoke: CZ", ncz, "url", box["url"], "bytes", len(json.dumps(box["payload"])))
    elif mode == "head":
        groups = [tests[i:i + 10] for i in range(0, len(tests), 10)]
        circs, plan, tot_cz = [], [], 0
        for g in [cal_a] + groups + [cal_b]:
            lq, ncz = pack(g)
            circs.append(lq); tot_cz += ncz
            plan.append([dict(j=t["j"], c=t["c"], w=t["w"], exact=float(t["u"] @ t["v"]), expect=t.get("expect")) for t in g])
        box = capture_payload(circs, shots=1000)
        json.dump(dict(plan=plan, pieces=pieces, q=q.tolist(), K=K.tolist(), **box), open("qpu_head_payload.json", "w"))
        print("head: circuits", len(circs), "total CZ", tot_cz, "url", box["url"], "bytes", len(json.dumps(box["payload"])))


def build_head_single(shots=600):
    q, K, pieces = model_case()
    tests = build_tests(q, K)
    u_ref = tests[0]["u"]
    cal_a, cal_b = calibration_tests(u_ref, +1.0), calibration_tests(u_ref, -1.0)
    groups = [tests[i:i + 10] for i in range(0, len(tests), 10)]
    payloads, plan, tot_cz = [], [], 0
    for g in [cal_a] + groups + [cal_b]:
        lq, ncz = pack(g)
        tot_cz += ncz
        box = capture_payload([lq], shots=shots)
        payloads.append(box["payload"])
        plan.append([dict(j=t["j"], c=t["c"], w=t["w"], exact=float(t["u"] @ t["v"]), expect=t.get("expect")) for t in g])
    json.dump(dict(plan=plan, pieces=pieces, q=q.tolist(), K=K.tolist(), payloads=payloads), open("qpu_head_single.json", "w"))
    print("single: tasks", len(payloads), "total CZ", tot_cz, "shots", shots)


if __name__ == "__main__" and sys.argv[1] == "head_single":
    build_head_single(600)
