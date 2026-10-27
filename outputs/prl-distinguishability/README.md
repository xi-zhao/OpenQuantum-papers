# Noisy-state distinguishability — revised manuscript

本目录的有效稿件是 `paper.tex` / `paper.pdf`。题目为 **Information spectrum and design dependence of noisy-state distinguishability**。

新版保留可证明的结论：任意 Pauli 分布的有限尺寸双侧界、固定 i.i.d. 噪声的阈值两侧指数分离、非退化临界信道族的高斯标度，以及退极化阈值附近由上下界夹出的首阶二次速率。它同时给出退相干下不同 2-design 的精确反例，以及稳定子态排除通用紧衰减指数的界。隐私放大指数按已有结果引用；不宣称发现物理冻结相变，不使用高斯替代模型作数值证据，不保证 PRL 接收或 hep-th 适配。

## 文件

- `paper.tex`：英文正文及内嵌参考文献。
- `paper.pdf`：由当前正文编译的短篇双栏 PDF。
- `figures/bounds.pdf`、`figures/dephasing.pdf`：论文矢量图；同名 PNG 为预览。
- `code/reproduce.py`：唯一有效的数值复现入口，无导入副作用。
- `data/haar_n8_raw.csv`：64 对独立 Haar 态，n=8，15 个噪声强度；保留每个样本值。
- `data/haar_n8_summary.csv`：均值、样本标准误及严格有限尺寸界。
- `data/critical_brackets.csv`：n=1024、65536 的信息谱上下界，非态模拟结果。
- `data/dephasing_exact.csv`：MUB 与 Haar 的解析退相干迹距离。
- `data/verification.json`：检查结果、随机种子、环境版本和数据摘要。
- `REVISION_NOTES.md`：审稿问题的逐项处理和仍存在的限制。
- `VALIDATION.md`：数值、排版、复现验证范围及已知构建提示。
- `REVIEW.md`：对初稿的审稿意见，链接指向历史版本。
- `archive/initial-draft/`：原稿、原始启发式脚本和审稿时证据的只读历史快照；不属于新版复现链。

## 复现

已使用现有 Python 3、NumPy、SciPy、Matplotlib、TeX Live 与 REVTeX 环境；没有安装新依赖。确切 Python 包版本写入 `data/verification.json`。

从本目录运行以下任一命令。脚本以自身路径定位数据，所以从其他目录以绝对路径调用也可。

只读运行独立检查，不生成文件：

```bash
python3 -I code/reproduce.py --check-only
```

利用已保存的原始样本重新生成全部表格与图：

```bash
python3 -I code/reproduce.py
```

使用固定种子重新计算全部 Haar 原始样本（会更新本任务的生成物）：

```bash
python3 -I code/reproduce.py --resimulate
```

编译论文，辅助文件保留在忽略的 `build/` 中：

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper.tex
```

成功后更新交付 PDF：

```bash
cp build/paper.pdf paper.pdf
```

数值脚本检查显式 Pauli 求和与逐比特通道实现、n=1/2 完整 MUB 的二阶矩与解析迹距离、相关两比特 Pauli 分布的精确二阶矩和双侧界、信息谱 CDF、Rényi 优化恒等式及完全退极化极限。样本均值落在界内仅是有限尺寸一致性检查，不能替代证明，也不是对每个样本的保证。

图 1(a) 的噪声扫描复用相同态对，不同噪声点因而相关；图中误差条只表示每个噪声点的样本标准误。图 1(b) 完全由严格界计算，不把大 n 结果伪装成精确态模拟。图 2 是解析值，无抽样误差。
