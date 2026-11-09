"""生成 QZ03 量子随机数线路的提交 JSON（官方 SDK 序列化，离线，不含 API key）。"""
import json, sys
import qpu_head as H
from lqcloud import QuantumCircuit as LQC

n = int(sys.argv[1]) if len(sys.argv) > 1 else 35
shots = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
qc = LQC(n, n)
for i in range(n):
    qc.h(i)
qc.barrier()
for i in range(n):
    qc.measure(i, i)
box = {}
class _Stop(Exception): pass
def fake(session, method, url, **kw):
    box.update(url=url, payload=kw.get("json")); raise _Stop()
H.cp.request_with_retry = fake
prov = H.LQCloudProvider(api_key="offline-placeholder", url="https://cloud.logicalqubit.com", interactive=False)
prov.get_backend_config = lambda name: H.QZ03
be = prov.get_backend("QZ03", verify=False); be._config_cache = H.QZ03
try:
    be.run(qc, shots=shots, initial_layout=list(range(n)), result_format="memory")
except _Stop:
    pass
p = box["payload"]
json.dump(p, open("qrng_payload.json", "w"), separators=(",", ":"))
c = p["command"]["circuit"]
print("url", box["url"], "| n_qubits", c["n_qubits"], "| shots", c["shots"], "| format", c["result_format"], "| instr", len(c["instructions"]), "| bytes", len(json.dumps(p, separators=(',', ':'))))
print(sorted({i["name"] for i in c["instructions"]}))
