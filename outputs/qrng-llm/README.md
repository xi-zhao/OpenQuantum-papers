# 量子真机随机数池驱动 LLM 解码 —— 论文与可复现代码

把云端门模型量子真机（逻辑比特云 QZ03，35 个物理比特、未纠错 NISQ）当作大语言模型解码阶段的物理随机数源，
端到端驱动真实预训练模型 stories15M 的采样，并完成随机性、安全、吞吐与成本的全链路评测。

## 目录
- `paper/main.pdf`、`paper/main.tex` —— 中文论文（XeLaTeX + Noto CJK 编译；本机没有 Noto CJK 时自动改用 macOS 自带的宋体、冬青黑体）。配图 `fig_*.pdf`。
- `code/` —— 全部脚本：
  - `qrng_payload.py` 构造 35 比特 QRNG 线路；真机原始比特见 `data/qz03_qrng_raw.json`。
  - `nist_sts.py` NIST SP 800-22 子集检验；`analyze_qrng.py` 偏置/最小熵 + Toeplitz 提取 + 检验。
  - `qsample.py` 把随机数池接入 stories15M 采样，对照 MT19937 / CSPRNG。
  - `attack.py` MT19937 种子劫持/预测攻击攻防演示。
  - `throughput.py` 吞吐/延迟/成本分析。`figures.py` 生成论文配图。
  - `model.py` stories15M 加载/前向/分词；`qfull.py`、`qpu_head.py` 为附录"全前向量子化"的实现。
- `data/` —— 真机原始比特与各步结果 JSON，以及提取后的比特 `extracted_bits.npy`。

## 复现
脚本在**当前目录**读写数据，所以先进入 `data/`，再调用 `../code/` 下的脚本；`figures.py` 会把配图直接写进 `paper/`。

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib   # 附录"全前向量子化"另需 qiskit qiskit-aer 与逻辑比特云 SDK lqcloud
cd data
# 1) 随机性分析（输入 qz03_qrng_raw.json）
python ../code/analyze_qrng.py
# 2) 采样对照、攻防、吞吐（前两个需要模型文件，见下）
python ../code/qsample.py && python ../code/attack.py && python ../code/throughput.py
# 3) 配图（写入 ../paper/）与论文
python ../code/figures.py
cd ../paper && xelatex main.tex && xelatex main.tex
```

## 模型文件（未随包提供）
`model.py` 读取 `stories15M.bin` 与 `tokenizer.bin`（Karpathy 的 llama2.c 格式，公开可得），
默认放在项目根目录的 `model/` 下，也可用环境变量 `QRNG_LLM_MODEL_DIR` 指向其他目录。
- `stories15M.bin`：https://huggingface.co/karpathy/tinyllamas/resolve/main/stories15M.bin
- `tokenizer.bin`：https://github.com/karpathy/llama2.c/raw/master/tokenizer.bin

## 真机平台
逻辑比特云 `cloud.logicalqubit.com`，SDK `lqcloud`。
**安全提醒：本包已移除任何 API 密钥。** 由于此前曾在会话中使用过你的 API key，
建议你到控制台\*\*轮换（重置）该 key\*\*。提交请求时把 key 放在环境变量或本地配置里，不要写进代码或数据文件。

## 诚实的限定
- 35 比特 QZ03 电路可被经典模拟，故此处随机性不具备 JPMorganChase–Quantinuum 式"可认证随机性"的密码学保证；
  随机性保证来自 NIST 统计检验与提取器的剩余哈希引理。
- 成本（约 ¥364/Mbit）远高于专用光子 QRNG；本方案价值在于证明通用云端量子机可胜任、且为无损即插即用。
