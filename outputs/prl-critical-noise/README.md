# 临界态首光子统计：研究记录与创新性否决

**当前判断（2026-10-29）：不作为 PRL 主选题。** 已完成可复现的解析与数值工作，但最主要机制及完整强监测生存函数与 Fedorov et al., JETP97,522(2003) 精确重合。本目录没有交付一个达到用户要求的 PRL 级新成果。

主要文件：

- [英文短篇研究记录](paper/research_note.pdf)：≤5页，明确标注创新限制；对应 [LaTeX](paper/research_note.tex)。
- [中文创新性审查](review/NOVELTY_AUDIT.md)：正式原刊 Eqs.69、78、84 与当前计算的精确变量对照，已读与未读文献分开列明。
- [独立审查记录](review/INDEPENDENT_REVIEW.md)：解析核验、Gaussian反例、经典随机过程表示与正式全文复核。
- [数值验证](data/verification.json)：截断、网格/步长、概率归一化、旧公式归一化检查。
- [第一幅图](figures/first_click.png) 与 [微观收敛图](figures/microscopic.png)。
- [图中数值表](data/figure_table.tsv)。
- [文献来源](sources/README.md)：下载全文仅本地保留，不纳入源码包或Git。

## 结果范围

有效模型 h=P²/2+Q⁴/4，jump=√d Q，偶基态初态。确定性无点击传播支持中位数∝d^−1、均值∝d^−3/4及分段矩标度。原始Dicke模型在N=16–512、固定d=.3/1/3下趋向软模预测；最大尺寸仍有可见有限尺寸修正，不能称精确收敛。微观保存的truncated_mean不是完整均值。

强监测边界层K(u)已在2003正式论文Eq84给出。普通Gaussian初态也有相同指数，且领先Wigner方程有经典随机加速加吸收表示。因此不得将其描述为新的量子临界普适类、本质量子等待时间异常或非经典见证。

只有旋转框架下工程化裸腔loss模型被考虑；结论不能用于保持真实基态不变的零温dressed bath。没有探测低效率、暗计数、实际硬件、完整Lindblad纯度动力学或微观矩统一尾界证明。

## 复现

已有Python3.13.2、NumPy2.2.6、SciPy1.16.3、Matplotlib3.10.7及TinyTeX/REVTeX。没有新增依赖。以下命令从仓库根目录运行。

从保存数据重建图和验证报告：

```bash
OPENBLAS_NUM_THREADS=1 python3 -I outputs/prl-critical-noise/code/analyze.py
```

重新计算最终偶基态矩：

```bash
OPENBLAS_NUM_THREADS=1 python3 -I outputs/prl-critical-noise/code/first_click.py --d 10 100 1000 10000 100000 --grid 9601 --dt-factor 0.002 --output soft_mode_moments.json
```

奇初态对照：

```bash
OPENBLAS_NUM_THREADS=1 python3 -I outputs/prl-critical-noise/code/first_click.py --d 10 100 1000 10000 100000 --grid 4801 --dt-factor 0.002 --state 1 --output soft_mode_odd.json
```

微观截断对照：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -I outputs/prl-critical-noise/code/dicke_first_click.py --N 256 512 --d 0.3 1 3 --cutoff 64 --output dicke_cutoff_check.json
```

重新生成带SHA-256清单的源码包：

```bash
python3 -I outputs/prl-critical-noise/code/package_source.py
```

编译短稿：

```bash
latexmk -pdf -cd -interaction=nonstopmode -halt-on-error -outdir=build outputs/prl-critical-noise/paper/research_note.tex
```

生成物位于paper/build；交付PDF是核验后复制到paper的副本。修改源文件后需要同步。

## 记录解释

- `soft_mode.json`：最初12个d点，grid1201、step factor .008。
- `soft_mode_refined.json`：grid2401、step factor .004。
- `soft_mode_converged.json`：文件名是历史命名，grid9601、step factor .002；不表示严格误差已为零。
- `soft_mode_moments.json`：同样细网格，补充发射通量与分数矩。
- `soft_mode_odd.json`：state=1的局部节点对照；JSON顶层constants仍指偶基态，不能用于奇态渐近系数。
- `dicke.json` / `dicke_large.json` / `dicke_cutoff_check.json`：N16–128 cutoff36；N256–512 cutoff48；同样大N cutoff64。

后续研究若复用本目录，必须保留已确认的先例边界，不能重新把K(u)、二次吸收节点或其长尾截断作为原创。
