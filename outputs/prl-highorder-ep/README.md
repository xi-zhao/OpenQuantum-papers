# PRA 短文：Liouvillian EP 的精确谱扇区与大尺寸标度

**稿件：** *Spectral sectors and finite-size scaling of fermionic Liouvillian exceptional points*。

- 阅读版：[`paper/manuscript.pdf`](paper/manuscript.pdf)，REVTeX/PRA 双栏，含图与参考文献共 4 页。
- 源码：[`paper/manuscript.tex`](paper/manuscript.tex)、[`paper/references.bib`](paper/references.bib)。作者、单位未填写，避免代填不明信息。
- 查新证据：[`review/prior_art.md`](review/prior_art.md)。正式论文元数据已核实；Xu–Yi 的逐式核对依据 arXiv v2，未逐页核对受限访问的 APS 出版 PDF。
- 验证数据：`data/verified_*.json`；图的数值表：`data/scaling_table.csv`、`data/splitting_table.csv`。
- 两幅图：`figs/splitting.pdf`、`figs/scaling.pdf`，均为矢量格式，另有 PNG 预览。

这是一篇可供作者审阅的 PRA 常规短文草稿。没有投稿、上传预印本或推送远端；完成写作和内部核验不等于外部同行评审或期刊录用。

## 核心论证

**1. Eq. (54) 的边界。** Xu–Yi 的 `z^(1/(d+1))` 指数仍然成立。其一阶特征多项式在固定非零谱参数下可以作为截断，但在分裂尺度 `μ~z^(1/(d+1))` 上，所有允许的高次项均为同阶。对 n=7,d=2，精确式为

`μ^8 − 504 z μ^5 + 27504 z² μ²`。

精确结果为两组三元环与一个二阶残余 Jordan 块；直接求一阶截断的根会得到一组三元环和五个零根。高次项首次出现的条件是 `n+1>=2(d+1)`；没有据此否定源文 n=4、d=2,3,4 的图。

**2. 全谱归约及宇称。** 对任意粒子数守恒的二次 Hamiltonian，全部等权 d 子集产生通道允许先旋转单粒子基，再按 ket-only/bra-only 冻结集合分块。消去费米符号后，每块是集体自旋 `J_-+z J_+^d/d!`，进一步化为不同尺寸的 `N`。精确谱隙须比较最大块的内部间距与相邻允许自旋块的 Perron 根差。全算符空间、物理偶宇称空间和总粒子数所见的间距不能混用。

**3. 固定 d 的半经典展开。** 定义

- `α=2d/(d+1)`，`β=α−1`；
- `a_d=[2^(−d+1)/(d−1)!]^(1/(d+1))`，`C_d=a_d/α`；
- `Ω_d=sqrt((d+1)/2)`。

大尺寸边缘阶梯为

`μ_l=z^(1/(d+1)){C_d n^α + a_d[1−(2l+1)Ω_d]n^β + O(n^(β−1))}`。

这是形式 WKB/大自旋展开；局部阶梯与前两正径向根的识别通过数值验证，尚未给出统一谱余项界。正文将其与精确恒等式分开陈述。

**4. 后选择的物理条件。** 原 recycling-only 生成元迹不增要求 `z<=n/binom(n,d)`。固定 d 的大 n 物理缩放可选 `z=ζ n^(1−d)`、`ζ<d!`；谱隙可趋有限，而存活概率的渐近指数仍随 n 广延衰减。稿件不主张传感优势或低成本后选择。

## 新意和先例的边界

| 内容 | 稿件中的定位 |
|---|---|
| 结构化 Jordan 扰动、多个环和全同阶多项式 | Ma–Edelman (1998)、Graefe 等 (2008) 已有 |
| 正双带矩阵的精确射线谱、非零单根 | McMillen (2009) 定理 |
| d=1 自旋精确解 | 源文 Appendix C 及既有结果 |
| d=2 领先谱端常数 C₂ | Shapiro–Tater (2022) 已证明；本文恢复该结果 |
| 本模型的完整冻结集合/自旋归约与宇称分辨公式 | 本稿明确给出，带矩阵恒等式验证 |
| 一般固定 d 的首个有限尺寸修正、局部间距及物理缩放 | 本稿半经典结果；未见直接近邻给出同一公式，不声称穷尽查新 |

## 验证范围

1. **精确代数**：6 个非平凡 `(n,d)` 组合的符号特征多项式、零根代数/几何重数；矩阵层面的 z 齐次相似变换；d=1 已知谱。
2. **Perron 根**：d=2,3,4，n=16 至 512 的 11 个尺寸及相邻 n−1,n−2，共 99 次 log-space 迭代。全部达到相对 Collatz 区间宽度 `<2e-12`。这是浮点收敛诊断，不是严格区间算术。
3. **多精度谱隙**：d=2,3,4；n=16,32,48,64；分别用 50 与 70 位精度独立解谱，间距差 `<1e-10`，同时核对正迭代 Perron 根。
4. **完整矩阵**：对 1≤n≤4 的所有 1≤d≤n，分别使用随机复 Hermitian 单粒子 Hamiltonian 验证基变换与符号规范后的矩阵恒等式，最大残差约 `2e-14`。另有 4,134 项整数 CAR 符号检查、自旋 intertwiner 至 m=7，以及 49 项主子矩阵恒等式。
5. **反例、解析极限与物理边界**：n=d=4 证明内部谱模式有时控制 gap；d=1 的全算符和偶宇称 gap 分别等于 `γ sqrt(z)` 与 `2γ sqrt(z)`。另对三种模型直接计算 `L†I`，核验迹不增区域内、边界上和边界外共 9 个点。

完整报告保留在 `data/verified_algebra.json`、`verified_perron.json`、`verified_precision.json`、`verified_blocks.json`。有审阅证据的结果见 `review/`，不依赖脆弱的高阶 EP 数值本征值逐根配对。

## 复现

以下命令从本目录执行。使用已有 Python 3.13.2、NumPy 2.2.6、SciPy 1.16.3、mpmath 1.3.0、SymPy 和 Matplotlib。未安装新依赖。环境版本进一步记录在 `data/verification_summary.json`。无需网络请求、模型账户或实际量子硬件。

```bash
python3 code/spectral_analysis.py algebra
```

```bash
python3 code/spectral_analysis.py perron
```

```bash
python3 code/spectral_analysis.py precision
```

```bash
python3 -I review/block_parity_review.py > data/verified_blocks.json
```

```bash
python3 code/make_figures.py
```

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/manuscript.tex
```

多精度计算比其他步骤慢，具体耗时记录在 JSON provenance 中。正迭代不收敛时会明确报错，不输出伪成功数据。全文均用实际数值，无拟合常数。图采用固定颜色、不同标记、文字标签及附带数值表，适合黑白打印。

## 探索记录隔离

旧的 `code/t*.py`、`perron2.py`、`figdata.py`、`wkb.py`、`data/perron_mu.json` 保留在原目录作为探索历史，**不作为稿件证据、不纳入正式源码包或本任务提交**。已知问题包括非正规矩阵双精度失真、未收敛的普通幂迭代、早期错误的 WKB 常数，以及 t19 中退化/多分支填充率方程。正式入口仅为 `spectral_analysis.py`、`make_figures.py` 和独立矩阵复核脚本。更正后的预检说明见 `PRECHECK_B.md`。
