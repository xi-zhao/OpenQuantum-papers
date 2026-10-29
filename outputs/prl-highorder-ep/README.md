# PRA 短文：费米 Liouvillian EP 附近的观测量弛豫

**稿件：** *Observable-dependent relaxation near a fermionic Liouvillian exceptional point*。

- [主文 PDF](paper/manuscript.pdf)：PRA/REVTeX 双栏，含两幅图和参考文献共 4 页。
- [补充材料 PDF](paper/supplement.pdf)：独立 4 页，保留完整分块、可观测量投影、谱多项式和半经典推导。补充页数不计入上述主文 4 页。
- [主文源码](paper/manuscript.tex)、[补充材料源码](paper/supplement.tex)、[参考文献](paper/references.bib)。作者、单位留空，供作者填写。
- [全文查新记录](review/prior_art.md)、[可观测量选择审阅](review/observable_selection_review.md)、[半经典边界审阅](review/semiclassical_review.md)。这些是内部推导检查，不是外部人类同行评审。
- 新图为 `figs/observable_dynamics.pdf`、`figs/rate_hierarchy.pdf`，均有 PNG 预览和配套 CSV。旧谱图保留在补充材料中。

本版按用户授权，将中心问题重构为：**同一生成元、同一个合法初态，为什么总粒子数和模式占据差看到不同弛豫率？** Eq. (54) 的截断问题已降为补充说明。完成写作和内部核验不代表期刊接受；未投稿、上传预印本或推送远端。

## 核心物理结果

在 Xu–Yi 的均匀 d 费米子 recycling-only pump 模型中，选择能量本征模式的

`D=n1−n2`，`rho0=p^N(1−p)^(n−N)[1+epsilon D]`。

对 `0<p<1`、`0<|epsilon|<1`，该态归一化、满秩正定、宇称为偶；它可以同时含奇数和偶数粒子。取 `p=.75, epsilon=.6`，两个可观测量都具有非零信号，不需制备奇宇称密度算符。

**精确选择规则：** 迹和总粒子数 N 只探测 Boolean 自旋 n/2 的块；D 只探测自旋 n/2−1 的一个固定副本。两个有限人口矩阵即可精确给出条件化的 N、D 和存活概率。改变 epsilon 不影响 N 与存活率的整条曲线，只线性改变 D。

令 `q=z^(1/(d+1))`，`kappa_M` 为尺寸 M 块在 z=1 的 Perron 根，`delta_n` 为最大块内部的实部间距。精确谱隙为

`Delta_all=gamma q min{delta_n,kappa_n−kappa_(n−1)}`，

`Delta_even=gamma q min{delta_n,kappa_n−kappa_(n−2)}`。

当 `n−2>=d`，上述有偏正初态保证

`Gamma_D=gamma q[kappa_n−kappa_(n−2)]`。

N 的指数取决于归一化后实际可见的内部根；只有相应重叠非零时，才可写成 `Gamma_N=gamma q delta_n`。复根带振荡，Jordan 根可带时间多项式。小尺寸 `n−2<d` 时 D 的长时形式为 `t^(n−2) exp(−gamma q kappa_n t)`。不能把 Gamma_D 无条件等同于整个偶宇称谱隙；`n=d=4,z=1` 是物理范围内的反例。

主图例子 `n=12,d=3,z=1/144` 的结果为：

| 量 | 数值，单位 gamma |
|---|---:|
| 全算符谱隙 Delta_all | 0.571963 |
| 偶宇称谱隙 / D 的指数 | 1.117616 |
| 总粒子数的指数 | 1.678630 |

这些量从同一个初态比较。偶态的一时刻偶观测量不含奇扇区；**不把这一点扩展为奇扇区永远不可观测**，费米双时 Green function 可以有奇中间插入。

## 读出混合与后选择代价

对 `O_eta=N/n+eta D/n`，生成元与初态不变，只改变读出权重。若 Gamma_D<Gamma_N，任何固定非零 eta 的渐近尾部由 Gamma_D 控制；eta=0 保留较快的 Gamma_N。

这是普通的两指数竞争在本模型中的定量应用，不是新的一般原理。主图理论振幅由左右本征矢投影得到，无拟合参数；另做尾部拟合仅用来交叉检查。

对主图例子，eta=.1 的估计交叉时间为 `gamma t≈10.54`。但存活概率在 `gamma t=5` 已降至 `1.18e−15`，在 `gamma t=8` 为 `2.02e−25`。图中直接显示这一代价，**不主张易于实验观察的长时精度优势**。

原生成元迹不增要求 `z<=n/binom(n,d)`。固定 d 的大 n 物理缩放可选 `z=zeta n^(1−d)`、`0<zeta<d!`；条件化弛豫指数可趋有限，存活概率的长时衰减指数仍随 n 广延。

## 半经典预测与查新边界

固定 d、局部低激发的形式展开为

`mu_l=q{C_d n^alpha+a_d[1−(2l+1)Omega_d]n^beta+O(n^(beta−1))}`，

其中 `alpha=2d/(d+1)`、`beta=(d−1)/(d+1)`、`Omega_d=sqrt((d+1)/2)`、`a_d=[2^(−d+1)/(d−1)!]^(1/(d+1))`、`C_d=a_d/alpha`。

在径向谱端识别和可见重叠的条件下，预测

`Delta_all : Gamma_D : Gamma_N -> 1 : 2 : sqrt[2(d+1)]`。

这不是任意 n/d 的精确定理。局部阶梯与前两正径向根的识别有数值支持，尚未证明统一全局谱余项界。d=1 给出已知精确极限 1:2:2。

| 内容 | 定位 |
|---|---|
| 观测量只探测特定扇区、不同观测量可有不同衰减率 | van Caspel–Gritsev, PRA 97, 052106 (2018) 等已有 |
| 分数指数、结构化 Jordan 扰动与多环 | 源文、Ma–Edelman (1998)、Graefe 等 (2008) 已有 |
| 正双带精确射线谱、非零单根 | McMillen (2009) 定理 |
| d=2 领先谱端常数 | Shapiro–Tater (2022) 已有，本文恢复 |
| 本模型 N/D 的精确扇区归属、同态条件动力学、读出混合的定量结果 | 本版主线，附代数和传播验证 |
| 完整冻结集合归约与宇称分辨 min-gap 公式 | 精确代数结果，附完整矩阵恒等式验证 |
| 一般固定 d 的首修正及观测量速率层级 | 形式半经典结果，保留证明边界；未声称穷尽查新 |
| Eq. (54) 的非一致截断 | 补充说明；不推翻源文指数或其 n=4,d=2,3,4 图 |

Xu–Yi 的正式出版元数据已核实；逐式核对依据 arXiv v2，未逐页核对受限访问的 APS PDF。完整查新记录保留证据位置与持久链接。

## 验证范围与异常记录

- **观测量选择**：n=2,...,7 的所有 d，共 27 项整数 intertwiner 检查；最大自旋与选中单态副本的残差均严格为零。
- **真实密度矩阵传播**：对 (n,d)=(4,2),(4,3),(5,2)，完整 CAR 超算符与标记模式降维传播的最大误差为 `1.34e−15`。初态奇分量和偶奇混合矩阵元为零；保持正性。
- **大一点的动力学实例**：(8,2),(12,3),(16,4) 三组。独立的两通道公式与标记人口过程的 N/D 差异小于 `3.4e−16`。长时谱振幅在 gamma t∈[6,9] 与传播信号相差至多 .21%；独立尾部拟合的速率误差至多 .23%。
- **如实记录失败**：n=4,d=3 的完整偶块双精度本征值在残余 Jordan 块产生 `5.53e−6` 伪分裂，不能据此核验精确 gap。数据保留偏差，将这一直接本征值比较标为不适用，依据整数选择规则、完整矩阵恒等式及实际传播验证；没有放宽误差阈值宣称通过。
- **既有谱证据**：6 项符号多项式；99 次 Perron log-space 迭代至 n=512，相对 Collatz 宽度 `<2e−12`；12 组 50/70 位谱隙核对至 n=64，保存结果差异 `<1e−10`；完整随机复 Hamiltonian 分块检查最大残差 `1.96e−14`；4,134 项整数 CAR 符号、49 项主子矩阵恒等式和 9 个迹不增边界检查。
- 浮点 Collatz 界是收敛诊断，不是严格区间算术。多精度结果转双精度后比较，不将零差解释为 50 位有效数字认证。

正式数据入口为 `data/verified_*.json`，汇总见 `data/verification_summary.json`。

## 复现

以下命令从本目录执行。使用已有 Python 3.13.2、NumPy 2.2.6、SciPy 1.16.3、mpmath 1.3.0、SymPy 1.14.0、Matplotlib 3.10.7 与 REVTeX 4.2。未新增依赖，无需网络、模型账户或量子硬件。

本版新增验证和作图：

```bash
python3 -W error::RuntimeWarning code/observable_relaxation.py
```

```bash
python3 code/make_observable_figures.py
```

已有谱证据的复现入口：

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

编译和打包：

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd paper/manuscript.tex paper/supplement.tex
```

```bash
python3 code/package_source.py
```

多精度计算较慢，耗时记录在 JSON 中。打包使用显式文件清单，逐文件 SHA-256 校验，不打包下载文献或探索脚本。

## 探索记录隔离

旧的 `code/t*.py`、`perron2.py`、`figdata.py`、`wkb.py`、`data/perron_mu.json` 留作探索历史，不作为稿件证据、不纳入正式源码包或本任务提交。已知问题包括非正规矩阵双精度失真、普通幂迭代未收敛、早期错误的 WKB 常数，以及 t19 的填充率分支失败。正式复现仅依赖上述入口。更正后的预检说明见 [PRECHECK_B.md](PRECHECK_B.md)。
