# Classical and quantum recovery in random codes

当前稿件：**A Common Recovery Transition for Classical and Quantum Information in Random Codes**，2026-10-28。

同一个随机量子编码器用于两项任务：读取均匀逻辑基消息，或恢复与参考系的纠缠。接收方知道编码器与固定 Pauli 信道，分别优化解码。正文给出任意 exact unitary 2-design 的有限尺寸信息谱夹逼；在固定 iid 噪声、0 < H < 1、正 surprisal 方差的正码率窗口中，两个最优性能具有相同的高斯极限与二阶码长，并给出编码器概率界。

文献定位：hashing 率、量子二阶可达性及量子信息谱可达侧已有直接先例（包括 Klesse 2007）。当前结果的重点是经典最优读取的匹配 converse，以及与量子恢复的共同操作性能。它不等于一般 Pauli 信道容量定理，不承诺高效解码，也不包含含噪编码线路或实际硬件实验。

## 当前文件

- `paper.tex` / `paper.pdf`：英文短篇正文与内嵌参考文献。
- `supplement.tex` / `supplement.pdf`：完整证明、Clifford 最优性、数值方法和前作对照。
- `code/recovery.py`：当前数值复现入口；导入不会执行模拟或写文件。
- `figures/recovery.pdf` / `figures/recovery.png`：论文矢量图与预览。
- `data/recovery_clifford_raw.csv`：4、8、12、16 比特各 64 个 Clifford 编码器，12 个噪声强度，共 3072 条记录。
- `data/recovery_encoders.json`：每个样本的种子及逆编码 symplectic matrix。
- `data/recovery_clifford_summary.csv`：每个尺寸、噪声强度的均值与标准误。
- `data/recovery_critical_brackets.csv`：大 n 的解析信息谱界；不是量子解码模拟。
- `data/recovery_rates.csv`：共同随机码率与已知最优经典容量的比较。
- `data/recovery_verification.json`：原始/对偶 SDP、显式恢复、枚举和 CDF 检查，含求解器尝试、版本与文件摘要。
- `VALIDATION.md`：本轮验证范围与结果。
- `REVIEW_ROUND11.md`：本轮独立模型审稿与处理记录。
- `LITERATURE_ROUND11.md`：直接相关原始文献、精确结果位置与重叠范围。
- `recovery-manuscript-source.zip`：当前正文、补充材料、源码、图表、数据和审查记录的交付包；包内有 SHA-256 清单。
- `REVISION_NOTES.md`：累积修订历史；较早结论仅代表相应轮次。

旧两态研究的 `code/reproduce.py`、不以 `recovery_` 开头的数据与旧图保留为历史记录。`archive/round10-memory/` 保存重写前正文与 PDF。`lean/` 仅覆盖旧两态稿的条件形式化，**没有形式化当前恢复定理**。

## 复现

使用现有 Python、NumPy、SciPy、Qiskit、CVXPY、CLARABEL、SCS、Matplotlib 以及 TeX Live / REVTeX；本轮未安装依赖。准确版本见验证 JSON。以下命令从本目录运行，脚本也支持绝对路径调用。

只运行数值检查，不改写数据或图：

```bash
python3 -I code/recovery.py --check-only
```

利用保存的样本生成汇总表、图和验证报告：

```bash
python3 -I code/recovery.py
```

按固定种子重新计算全部 Clifford 样本，并更新生成物：

```bash
python3 -I code/recovery.py --simulate
```

编译正文：

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper.tex
```

编译补充材料：

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build supplement.tex
```

编译产物保存在 `build/`，目录根的同名 PDF 是交付副本。若自行修改源文件，请在检查后同步副本。

Clifford 数值对有限尺寸模型精确求和，仅有浮点误差；它不以独立随机向量代替同一编码器的相关错误像。噪声扫描复用同一编码器，因此不同噪声点相关；误差条是编码器样本均值的标准误。小系统数值是独立一致性检查，不能替代渐近定理的证明，也不能单独证明 PRL 重要性。
