# 科学审稿意见

审查对象：初稿 `paper.tex` / `paper.pdf`，2026-10-27 版本，4 页。该历史版本现保存在 `archive/initial-draft/`；本报告的论文与旧脚本链接均指向该快照。当前修订稿的处理情况见 [REVISION_NOTES.md](REVISION_NOTES.md)。

**结论：当前稿不适合按现有主张投稿 PRL，需要实质重构。** 主定理的有限尺寸双侧界正确，有可保留的数学内容；但临界点结论存在解析反例，“真实衰减指数与冻结相变”的论证不成立，相关指数和临界结构也有重要的未引用先例。当前内容自然属于 quant-ph，尚无实质 hep-th 模型或结论。

这一结论纠正了上一轮撰稿时过强的自评。数值落在严格界内不能证明界的紧性，未经证明的替代模型也不能验证原问题的渐近指数。本次未修改稿件或原始脚本。

## 可保留的结果

1. **Theorem 1 正确。** 设两个纯态独立抽自精确 state 2-design，Pauli 分布可任意相关。下界的支持子空间投影测量和上界的二阶矩计算均成立。原文的精确恒等式为
   \[
   \mathbb E\operatorname{tr}\mathcal N_j(\psi_1-\psi_2)^2
   =\frac{2d}{d+1}\sum_{P\in\mathcal C_j}p_P^2
   -\frac{2p(\mathcal C_j)^2}{d(d+1)}.
   \]
   因而上下界及 \(1/\sqrt2\) 常数可以保留。下界实际上只需第二态构成 1-design，第一态可为任意纯态。
2. **对固定 i.i.d. Pauli 信道，\(H(p)<1\) 时 \(\bar T\to1\)、\(H(p)>1\) 时 \(\bar T\to0\)，且误差指数小，成立。** 这不包括 \(H=1\) 的所有退化情形。
3. **非退化临界点附近的高斯轮廓可成立。** 对规则信道族，若 \(\sigma\) 不趋于零且 Berry–Esseen/Lindeberg 条件一致成立，选择 \(a=K\log_2n\) 即可由定理得到 \(\Phi(x)+O(\log n/\sqrt n)\)。原稿选择 \(a=n^{1/4}\) 能给极限，但不能直接得到图注声称的这个误差阶。
4. **类型法给出的 \(E(p)\) 是有效的迹距离上界指数。** 若定义真实衰减率 \(R_n=-n^{-1}\log_2\bar T_n\)，正确方向是 \(\liminf R_n\ge E(p)\)，不是已证明 \(R_n\to E(p)\)。在非退化去极化临界点附近，配合原文下界可以夹出首阶二次行为；一般参数下仍未证明匹配。

## 主要问题

### 1. 临界点“所有 2-design 都有 T̄→1/2”是假的

**位置：** [paper.tex:60](archive/initial-draft/paper.tex:60)、[paper.tex:65](archive/initial-draft/paper.tex:65)，并影响摘要和第 55 行的普适性表述。

取单比特完全 Z 退相干信道
\[
p=(1/2,0,0,1/2),\qquad H(p)=1,\quad\sigma(p)=0.
\]
在 \(d=2^n\) 维，取包含计算基的完整 \(d+1\) 组互无偏基（MUB），对全部 \(d(d+1)\) 个态均匀抽样。这是精确纯态 2-design。完全退相干后：计算基态保持本身；其余各基的态全部变为 \(I/d\)。直接平均得到
\[
\bar T_{\rm MUB}
=\frac{(d-1)(2d+1)}{d(d+1)^2}\longrightarrow0.
\]
同一信道的 Haar 输入则满足
\[
\bar T_{\rm Haar}=\frac{d-1}{2d-1}\longrightarrow\frac12.
\]
所以这是**相同信道、相同精确 2-design 条件、不同临界极限**的反例，不只是公式除以零的记号问题。单比特六态枚举给 \(5/18\)，与 MUB 公式一致；\(n=8\) 时两者分别约为 0.007737 和 0.499022。

“along any sequence”还需要加强。仅逐个 n 要求 \(p_i>0\)、\(\sigma>0\) 不够：可取 \(\varepsilon_n=n^{-6}\)，令
\[
p_n=(1/2-\varepsilon_n+\delta_n,\varepsilon_n,\varepsilon_n,1/2-\varepsilon_n-\delta_n),
\]
选择 \(\delta_n>0\) 使 \(H(p_n)=1\)。此时 \(\delta_n=O(\sqrt{\varepsilon_n\log(1/\varepsilon_n)})\)，\(n\|p_n-p_0\|_1\to0\)。乘积信道因而趋近完全退相干信道，沿 MUB ensemble 仍有 \(\bar T\to0\)，但每个 n 的 x 都为零且方差为正。

**必要修改：** 明确限制到非退化、满足一致中心极限定理条件的信道族；单独处理零方差临界信道。原文的去极化临界点方差非零，不受上述反例直接否定。

### 2. 2-design 不能保证稿件声称的同一紧衰减指数

**位置：** [paper.tex:55](archive/initial-draft/paper.tex:55)、[paper.tex:87](archive/initial-draft/paper.tex:87)、[paper.tex:101](archive/initial-draft/paper.tex:101)。

令 \(\lambda=1-r\)，均匀随机稳定子态满足
\[
\rho=\frac1d\left[I+\sum_{P\ne I}s_PP\right],\qquad
s_P\in\{0,\pm1\},\qquad\mathbb E|s_P|=\frac1{d+1}.
\]
局域去极化将每项乘以 \(\lambda^{w(P)}\)。由迹范数三角不等式，对两个独立随机稳定子态
\[
\bar T_n
\le\frac1{d+1}\sum_{P\ne I}\lambda^{w(P)}
=\frac{(1+3\lambda)^n-1}{2^n+1}.
\]
因此在 \(r=0.9\) 时
\[
\liminf_{n\to\infty}-\frac1n\log_2\bar T_n
\ge1-\log_2(1.3)=0.621488\ldots,
\]
而稿件公式给
\[
E(0.9)=\tfrac12[1-\log_2(1.03)]=0.478678\ldots.
\]
两者严格不相等。随机稳定子态是精确 state 2-design，所以第 55 行“All statements below ... apply verbatim”不能成立。

**必要修改：** 区分对所有 2-design 成立的界、在正则条件下的阈值标度、以及仅针对 Haar ensemble 提出的指数猜想。此反例不否定 Theorem 1，也不否定 E 作为上界指数；它否定的是普遍紧性。

### 3. 高斯替代模型不能支撑原模型的紧指数或“冻结相变”

**位置：** [paper.tex:95](archive/initial-draft/paper.tex:95)、[paper.tex:101](archive/initial-draft/paper.tex:101)，影响摘要、图 2 和讨论。

\(P|\psi\rangle\) 之间的小两两重叠，并不意味着整个 Pauli 轨道可以当作独立高斯列。轨道满足强全局约束，例如
\[
\frac1{d^2}\sum_PP\psi P=I/d
\]
对每个输入态都精确成立。独立高斯替代破坏了这个恒等式。

- **完全退极化 r=1：** 两个输出严格均为 \(I/d\)，故 \(\bar T_n=0\) 对每个 n 成立；稿件高斯模型给出正值，其预测指数为 1/2。真实衰减率在扩展实数意义上为无穷大。这直接排除了“所有 p 都有真实指数 E”的字面猜想。
- **非端点 r=0.9、n=10：** 真实模型由精确二阶矩得到
  \[
  \bar T\le\frac12\sqrt{\frac{2}{d+1}\big[(1+3\lambda^2)^n-1\big]}
  =0.01295238\ldots.
  \]
  我用自适应积分复核稿件 Eq. (7)，得到 0.02174021，仍明显超过真实模型的严格上界。这不是积分精度能修复的偏差。

原稿已经承认某些参数下违反上界，却随即称“所以只用它来取指数”。这需要一个新的论证：例如证明替代误差对指数不产生影响。当前稿没有该论证。n=8、10 两个尺寸的有限比值也不能推出渐近 \(O(1)\) 前因子。

**必要修改：** 将模型准确标为启发式；不能把其大 n 外推当作原问题的独立数值验证。没有匹配下界前，\(r^\star\) 只能称为**该上界的优化分支切换点**，不能宣称真实迹距离存在冻结相变。Haar 情形的内点指数猜想并未被以上端点反例彻底排除，但远未得到证明。

### 4. 创新性定位遗漏了直接相关的既有结果

**位置：** [paper.tex:26](archive/initial-draft/paper.tex:26)、[paper.tex:28](archive/initial-draft/paper.tex:28)、[paper.tex:76](archive/initial-draft/paper.tex:76)、[paper.tex:87](archive/initial-draft/paper.tex:87)。

#### 4.1 “H>1 时排除代数衰减区”可由 Liu 自己使用的工具直接推出

稿件复述 Liu 正文对中间相的描述基本准确，但其定理实际上更强。[Liu Theorem 3, Eq. (10)](https://arxiv.org/html/2606.02122v1#Thmtheorem3) 和 [Lemma 1, Eq. (31)](https://arxiv.org/html/2606.02122v1#Thmlemma1) 给
\[
\bar T_n\le2^{-H_2^\delta(\mathcal N^{\otimes n})/2}+4\delta,
\qquad
H_2^\delta\ge nh-C\sqrt{n\log_2(2/\delta^2)}.
\]
对固定 \(h=H(\mathcal N)>0\)，取 \(\delta=2^{-\kappa n}\) 且 \(\kappa>0\) 充分小，则充分大 n 时
\[
H_2^\delta\ge nh/2,
\qquad\bar T_n\le2^{-nh/4}+4\,2^{-\kappa n}.
\]
原始 [AEP Theorem 9](https://arxiv.org/html/0811.1221v3#Thmlemma9) 的条件 \(n\ge(8/5)\log_2(2/\delta^2)\) 仍允许 \(\kappa<5/16\)。Pauli 信道对应 \(h=H(p)-1\)。因此固定信道严格在阈值以上时，指数衰减是已有界的直接推论。

纠正原文中间相解释仍有价值，但不能将其包装成需要新工具才解决的理论缺口。这一推论不涵盖 h=0，也不自动消除 approximate-design 误差。

#### 4.2 变分指数与临界点有既有隐私放大结果

原稿的公式严格等价于
\[
E(p)=\max_{0\le s\le1}\frac{s}{1+s}\,[H_{1+s}(p)-1].
\]
该恒等式可将 \([H(t)-1]_+/2\) 写成 \(\max_{0\le u\le1/2}u[H(t)-1]\)，交换凸–凹极值后对 t 最小化，再令 \(s=u/(1-u)\) 得到。

这正是 [Hayashi, Eq. (18)](https://arxiv.org/html/1010.1358v5#S3.E18) 的 universal₂ 隐私放大指数在 R=1 的形式。其 [Eqs. (26)–(27)](https://arxiv.org/html/1010.1358v5#S3.E26) 已包含临界率
\[
R_{\rm crit}=H\!\left(p_i^2/\sum_jp_j^2\right),
\]
对应稿件的 escort 条件和碰撞熵分支。这里的“冻结点”在数学上是优化参数到达 s=1 端点，不能仅凭 REM 类比主张它是新的物理相变。

[Dupuis, Theorem 8](https://arxiv.org/html/2105.05342v2#Thmtheorem8) 还可直接用于 Pauli 轨道的随机输出信道，得同一指数上界。具体地，令
\[
\mathcal R^\psi(X)=\sum_PX_{PP}P\psi P,\qquad
\mathcal U(X)=\operatorname{tr}(X)I/d.
\]
2-design 二阶矩保证它是该文定义的 \(\lambda\le1\) randomizing family（含旁系）；该定理于是给
\[
\mathbb E_\psi\|\mathcal N_p(\psi)-I/d\|_1
\le2^{2/\alpha-1}
2^{-\frac{\alpha-1}{\alpha}[H_\alpha(p)-\log_2d]},\quad1<\alpha\le2.
\]
再对两态用三角不等式并优化 \(\alpha\)，即得到 E。上述既有结果也没有自动证明当前 pairwise observable 的真实指数等于 E。

#### 4.3 仍可谨慎保留的创新候选

[Sommers 等 §III–IV](https://arxiv.org/html/2510.07396v2#S4) 已讨论 hashing 阈值与 \(n^{-1/2}\) 转变窗宽：§IV 的“This picture of well-defined bands”段落明确给出窗宽，Fig. 3(b) 展示 coherent information 的标度塌缩，主文 Eq. (6) 给出其分段渐近式。但该观测量不是本文的同一两态迹距离。[Shen–Gao–Cheng](https://arxiv.org/html/2202.11590v1#Thmtheo10) 已有迹距离信息谱截断及隐私放大的二阶高斯展开，但也不能不经证明就转为本文相关 Pauli 轨道的精确 pairwise 高斯轮廓。

因此，**任意相关 Pauli 分布下的简洁双边界，以及正则条件下精确 pairwise 高斯轮廓，仍可能有独立价值。** 本次没有找到完整覆盖它们的现成定理；这不等于证明了新颖性，更不等于建立 PRL 所需的重要性。

### 5. 数值证据的误差与复现说明需要重做

**位置：** [code/theory2.py:15](archive/initial-draft/code/theory2.py:15)、[code/theory2.py:27](archive/initial-draft/code/theory2.py:27)、[code/figs.py:29](archive/initial-draft/code/figs.py:29)、[code/extra.py:15](archive/initial-draft/code/extra.py:15)、[paper.tex:91](archive/initial-draft/paper.tex:91)。

- 我把 `sim.py` 的逐比特通道实现与显式全部 Pauli 串求和在 n=1、2、3 比较，最大偏差分别为约 \(5.6\times10^{-17}\)、\(8.3\times10^{-17}\)、\(1.7\times10^{-16}\)。这部分实现没有发现错误。
- n=6、8、10 的真实密度矩阵样本只能检验这些尺寸。n=200、800 的图中点来自高斯替代模型，不能称为原信道的精确大尺寸数值证据；图 1(b) 也只画了替代模型，没有重画精确模拟样本。
- `theory2.py` 用均匀 log 网格等权求和，包含两个端点，没有积分收敛检查。在 n=640、r=0.5：

  | 积分方法 | \(-\log_2 T_{DE}/n\) |
  |---|---:|
  | 原算法 6001 点 | 0.121277330 |
  | 原算法 12001 点 | 0.121422629 |
  | 原算法 24001 点 | 0.121495813 |
  | 独立自适应积分 | 0.121569267 |

  图线所用网格产生的偏差约 \(2.9\times10^{-4}\)，已经和文中强调的 \(2\times10^{-4}\) 指数拟合差同量级。拟合差不能充当模型误差、积分误差和有限尺寸外推误差的合并误差条。
- 用自适应积分重做 n=100、200、400、800 的同一外推模型，r=0.35 的截距为 0.02162928，对照 E=0.02140518；r=0.5 为 0.11769933 对照 0.11758880；r=0.7 为 0.32759280 对照 0.32758575。这说明替代模型的拟合在部分点接近 E，但不验证原模型。
- `extra.py` 最后一行从 `scipy.misc` 导入已移除的 `derivative`；当前 SciPy 1.16.3 下无法成功结束。此前用管道或丢弃 stderr 运行，掩盖了失败状态。
- 脚本导入有副作用，部分会重新运行模拟或覆盖数据；随机 simplex 优化没有固定种子。尚无一次完整、退出状态受检查的复现入口。

### 6. 若干物理解释超出了证明，且有一句方向相反

**位置：** [paper.tex:55](archive/initial-draft/paper.tex:55)、[paper.tex:68](archive/initial-draft/paper.tex:68)、[paper.tex:103](archive/initial-draft/paper.tex:103)、[paper.tex:105](archive/initial-draft/paper.tex:105)。

- \(s<n\) 表示单个错误的概率大于 \(1/d\)，不是“最可能的 d 个错误”的同义语。该集合最多有 d 个元素，但可能为空；top-d 集合仍有 d 个元素。
- Theorem 1 是两个**不同阈值**的 CDF 夹逼，不是非渐近 identity。要近似 \(\Pr[s<n]\)，还须控制区间 \([n-a,n+a]\) 内的概率质量。
- 支持子空间投影只能检测一组错误像，不能识别每一个具体 Pauli 错误。小维数不自动证明可纠错，也不能证明未知量子态可恢复。
- 两个已知状态的二元判别与保持任意叠加态的量子传输不同；不能把 hashing 相干信息写成一般量子容量的充要条件。
- “Below the hashing point no noise-robustness argument based on scrambling can apply”与本稿 \(H<1\Rightarrow\bar T\to1\) 的结论方向相反。
- “matching constants on both sides”也不准确：两侧的指数误差和阈值移位不同。

## 期刊与主题判断

4 页满足用户的长度约束，但长度并不补足科学依据。当前尚不足以支持 PRL：两项被强调的主结果存在错误或未证明部分，剩余指数结构与标准理论有明确重合。比较稳妥的研究核心应是：**Pauli 噪声下随机态两两可区分度的有限尺寸信息谱刻画**，再说明它比现有 soft covering / decoupling 结果多提供了什么。

当前没有量子场论、引力、全息或黑洞模型；Hayden–Preskill 引用和 scrambling 术语不足以使论文成为 hep-th 结果。若最初 hep-th 是硬约束，这次选题尚未满足；若量子计算方向可以独立满足要求，quant-ph 是自然定位。PRL 本身接收量子信息论文，无须附加 hep-th 包装。

## 修订优先顺序

1. 保留 Theorem 1；修正 Corollary 的方差与序列正则性假设。
2. 将所有 2-design 的定理与 Haar 特有猜想分开；删除普遍紧指数与已发现物理冻结相变的断言。
3. 补引 AEP、隐私放大和 soft covering 的直接先例，重写创新性说明。
4. 用严格界、反例和受控的精确模拟呈现证据；高斯替代只作为明确标注的启发式，若无作用则删去。
5. 在以上完成后再评估是否有足以支持 PRL 的独立物理贡献。

这些是审稿建议，本次未执行修订。

## 核查范围与来源

核查了 LaTeX 原文、4 页 PDF 中的公式与图注、模拟和绘图脚本；独立复算了有限尺寸通道实现、解析反例、严格二阶矩界、Rényi 变分恒等式、积分收敛性和部分拟合。没有把此前小尺寸模拟通过当作新颖性或渐近正确性的证明，也没有重新运行所有耗时模拟。这些核查支持本报告中的具体审稿意见，不构成对其余未证明猜想的证明。

主要原始文献：

- [Liu et al., Can scrambling protect quantum state distinguishability under noise?](https://arxiv.org/html/2606.02122v1)
- [Tomamichel, Colbeck, Renner, A Fully Quantum Asymptotic Equipartition Property](https://arxiv.org/html/0811.1221v3)
- [Hayashi, Tight exponential analysis of universally composable privacy amplification and its applications](https://arxiv.org/html/1010.1358v5)
- [Dupuis, Privacy amplification and decoupling without smoothing](https://arxiv.org/html/2105.05342v2)
- [Shen, Gao, Cheng, Optimal Second-Order Rates for Quantum Soft Covering and Privacy Amplification](https://arxiv.org/html/2202.11590v1)
- [Sommers et al., Spectral properties and coding transitions of Haar-random quantum codes](https://arxiv.org/html/2510.07396v2)
- [PRL 官方范围及标准](https://journals.aps.org/prl/about)

审查版本指纹：

- `paper.tex` SHA-256: `b0a831dbd01c0a8e8f0fce06f3ec988f053c29757c94358bf2b00a6bc47200a2`
- `paper.pdf` SHA-256: `c6878c76af526318db1e60b2c24f580aeb1f2d7ee8652caa7da049ae957f3b20`
