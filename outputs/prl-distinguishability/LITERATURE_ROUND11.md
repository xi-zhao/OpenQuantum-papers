# 第十一轮：直接前作与贡献边界

核查日期：2026-10-28。对象为共同随机编码器的最优经典基消息成功率 P_c* 与最优纠缠恢复保真度 F_q*。本记录依据原始论文全文核查，不把检索摘要作为定理证据，也不把“未找到相同定理”写成不存在先例。

## 可辨认的增量

对抽样编码器之前固定的任意 Pauli 分布，当前稿件给出 exact unitary 2-design 系综上 P_c* 的信息谱 converse。结合量子可达界及逐编码器的 F_q* ≤ P_c*，确定正码率临界窗口内两个最优任务的共同 Gaussian profile，并控制同一个编码器样本的偏差概率。

Hashing 率、量子二阶可达性、2-design 随机纠错和高概率量子可达性均有先例。尤其是量子 surprisal-CDF 可达侧本身可由 Klesse 的旧结果推出；正文与补充材料已明确说明。

## 最直接的量子可达侧前作

R. Klesse, *Approximate quantum error correction, random codes, and quantum channel capacity*, Phys. Rev. A **75**, 062315 (2007)。[arXiv v2](https://arxiv.org/abs/quant-ph/0701102v2)，以下公式号对应该版本。

- Eq. (38) 对 trace-decreasing CP operation 给平均最优码保真度下界：tr N(π) − √(K L) ‖N(π)‖₂，其中 L 是所选 Kraus 表示的算子数；保真度未按截断成功概率重新归一化。
- Eq. (35) 保证删除部分 Kraus 项后，最优保真度不超过原信道。
- 应用于截断 Pauli 家族 G，令 q=p(G)、L=|G|，得到 E F_q* ≥ q[1−√(M L/d)]。取 surprisal < n−k−a，得到 E F_q* ≥ F(n−k−a)−2^(−a/2)。不存在 q² 前因子损失。
- 取 a=2 log₂ n 已足以推出当前正码率窗口的 Gaussian 量子可达侧，误差 O(log n/√n)。
- 原文陈述 Haar 随机码；Eq. (26) 只依赖码投影二阶矩，所以 exact unitary 2-design 延伸直接成立。不能写成原文已显式陈述 2-design 定理；Sommers v3 Sec. VI.1 也讨论这一延伸。

- Klesse 的逐码 Eq. (36) 还与平均估计合用给出 E(q−F_q*)₊ ≤ q√(ML/d) ≤ 2^(−a/2)。Markov 不等式即可推出量子侧高概率可达性；取 a=4 log₂ n 与 η=n^(−1/2)，仍保持 O(log n/√n) 的 Gaussian 精度和多项式消失的例外概率。

当前显式 square-root recovery 证明将罚项改进为 8·2^(−a)。主要增量是匹配经典 converse 和共同操作曲线；量子侧 Gaussian、2-design 与高概率不能分别计算为独立创新。

## 其他直接对照

| 前作与原始来源 | 核实的结论位置 | 与本稿的关系 |
|---|---|---|
| M. Tomamichel, M. Berta, J. M. Renes, *Quantum coding with finite resources*, Nat. Commun. **7**, 11419 (2016), [正式版](https://www.nature.com/articles/ncomms11419) | 正式版 Theorem 1、Eq. (6)：dephasing 的匹配 Gaussian/三阶展开；Theorem 5、Eqs. (15),(37)：一般量子二阶可达界；Eqs. (32)–(34)：一次性界及渐近化 | 误差采用纠缠 infidelity，直接相关。优化编码，不是指定任意 2-design 系综的两个最优任务共同曲线。正式版编号与 arXiv 的 Result 1/5 不同。 |
| A. S. Darmawan, Y. Nakata, S. Tamiya, H. Yamasaki, *Low-depth random Clifford circuits for quantum coding against Pauli noise using a tensor-network decoder*, Phys. Rev. Research **6**, 023055 (2024), [DOI](https://doi.org/10.1103/PhysRevResearch.6.023055) | 正式 Appendix D：D6 定义最大纠缠输入的 trace-distance 恢复误差；D10–D11 为 2-design decoupling；D20–D24 是二阶 smooth-entropy 可达估计及块长条件 | 单侧可达性；没有匹配 P_c* converse。正式 PDF 第 12–13 页包含此分析，arXiv 2212.05071v1 的同名附录不能代替。原文 D4/D18–D19 的 V 按其定义引用，不未经换算等同于本稿中心化 surprisal 方差。 |
| G. M. Sommers et al., *Spectral properties and coding transitions of Haar-random quantum codes*, PRX Quantum **7**, 020328 (2026), [v3 全文](https://arxiv.org/html/2510.07396v3) | Eqs. (6)–(8)：coherent information；Sec. IV、Fig. 3(b)：iid depolarizing 的 n^(−1/2) 展宽；补充 S49–S50：谱矩/熵；Sec. VI.1：既有 Hamming 可达界 | 随机码 hashing 转变及平方根窗口已有明确先例，主要讨论固定逻辑比特数的零码率情形。没有本稿最优 entanglement fidelity 与经典基消息成功率的共同有限尺寸夹逼。 |
| F. Ma, X. Tan, J. Wright, *Haar random codes attain the quantum Hamming bound, approximately*, [arXiv:2510.07158v1](https://arxiv.org/html/2510.07158v1) (2025) | Theorems 1.10–1.11、Definitions 1.2、1.5：固定 Hilbert–Schmidt 正交 unitary-error 家族的高概率近似纠错，半 diamond-norm disturbance | 量子误差指标及覆盖的相干 Kraus 组合强于当前量子可达侧，但限 Haar 且为可达结果；没有一般匹配 converse 或共同 Gaussian pair。“能纠正约 d/M 个 Pauli 错误”不能作为区别。 |
| M. Tomamichel, *A strong converse for stabilizer codes over Pauli channels via the blowing-up lemma*, [arXiv:2607.23450v1](https://arxiv.org/html/2607.23450v1) (2026) | Theorem 6.5、Eq. (58)：完整 stabilizer eigenspace 编码、k ≥ I_c(V)+γn 时最优 fidelity 指数衰减；Cor. 6.7、Eq. (61)：memoryless stabilizer-class strong converse | 阈值是码输入相关 coherent information / 正则化 Q_stab，不是普遍单字母 1−H(p)。Prop. 3.1、Eq. (23) 的 classical guessing 是 syndrome 给定下猜逻辑 Pauli 错误类，不是猜发出的逻辑基消息。 |
| Y. Nakata, E. Wakakuwa, H. Yamasaki, Phys. Rev. A **104**, 012408 (2021), [arXiv:2011.00668v2](https://arxiv.org/abs/2011.00668v2) | Theorem 1、Eqs. (16)–(18)：混合经典—量子消息的一次性恢复及资源区域，消息相关随机化编码 | 同时研究经典和量子传输本身已有先例。其任务不是同一个量子码的基消息读取与纠缠恢复最优值共同曲线。 |

## Liu 等：测量后互信息与优化顺序

G. Liu et al., *Can scrambling protect quantum state distinguishability under noise?*, [arXiv:2606.02122](https://arxiv.org/html/2606.02122) (2026)，Sec. IV.3、Appendix C，还讨论测量后的互信息、accessible information 和样本需求，不能把它只描述为二态问题。

其 state-ensemble 结论以整个 state-2-design 为字母表；本稿对每个固定且接收方已知的 U，字母表仅为 M 个正交列，之后优化这一条件字母表的 POVM，再对 U 平均。固定 U 的列集合本身不是 state 2-design。这里区别是字母表与优化顺序 E_U max_{POVM(U)}，而不只是消息数量。

## 判读范围

本轮核查未找到覆盖完整共同曲线的定理；这不证明所有文献都没有先例，也不能单凭定理成立证明 PRL 重要性。全文主张以当前 paper.tex 和 supplement.tex 为准。独立成稿审查的结论另记于 REVIEW_ROUND11.md。
