# PRA 稿件文献定位与引用建议

调研日期：2026-10-29。对象：*Finite-size decay-mode hybridization near a mean-field fold in an open Rydberg ensemble*，当前 [正文](paper/main.tex) 与 [PDF](paper/main.pdf)。本报告提供文献建议，不改动已交付稿件或其数值证据。

## 判断先行

目前的 12 条引用确实不够，问题主要在覆盖面：漏掉了**同模型的早期研究和精确稳态解、有限非零衰减率处的谱双重态、相反大残数相消，以及谱参数依赖消元的直接方法先例**。这些缺口会影响审稿人对我们是否了解现有工作的判断，也会改变哪些结果能够承担创新性。

建议优先新增下文 15 篇；现有 12 篇基本保留，核心参考文献由此变为 27 篇。时域记忆方程若继续保留，再加 Nakajima–Zwanzig 两篇历史出处；Kessler、Dutta、Macieszczak 2021 可按正文讨论深度加入。**27–32 篇是论证覆盖后的结果，不是需要凑足的指标。** [建议 BibTeX](paper/references_suggested_2026-10-29.bib) 包含 15 篇优先新增和 7 篇条件新增，未被接入正文。

本次核查未找到完整覆盖以下组合结果的论文：独立局域衰减的全连接 Rydberg 模型、完整置换不变算符空间、左 fold 单稳侧的有限率实模避交叉、物理乘积态扰动下高分支响应交换、以及该谱对在相干消元中的双自洽根结构。这是有边界的检索结论，不能写成“首次证明”或“已穷尽所有先例”。

## 最接近的模型与物理结果：优先新增 8 篇

### 1. Lee–Häffner–Cross 2012：模型的直接源头

**T. E. Lee, H. Häffner, and M. C. Cross, “Collective quantum jumps of Rydberg atoms,” Phys. Rev. Lett. 108, 023602 (2012).** [DOI](https://doi.org/10.1103/PhysRevLett.108.023602) · [arXiv](https://arxiv.org/abs/1109.6506) · [本地全文](src_dl/literature_2026-10-29_model/lee2012/paper.pdf)。

- **核对**：全文 Eqs. (1)–(2)、(6)–(7)，Figs. 3–5。激光驱动、全连接 Rydberg 相互作用、逐原子独立自发衰减；已使用 O(N³) 置换约化，研究平均场双稳态和有限尺寸量子切换。
- **与本文关系**：是同一有限 N 模型，差别仅为相互作用归一化；精确参数换算见后文。
- **落点**：Introduction 第一段；Model 的 Hamiltonian 与耗散项之后。Carr 2013 继续承担实验背景，不能代替这篇理论模型来源。
- **边界**：标题中的“collective quantum jumps”指多体协同切换，不意味着 Lindblad 跳跃是集体 J₋。未发现它计算本文所研究的有限率激发模间距。

### 2. Lee–Chan 2013：同模型平均场、Jacobian 与涨落分析

**T. E. Lee and C.-K. Chan, “Dissipative transverse-field Ising model: Steady-state correlations and spin squeezing,” Phys. Rev. A 88, 063811 (2013).** [DOI](https://doi.org/10.1103/PhysRevA.88.063811) · [arXiv](https://arxiv.org/abs/1309.3255) · [本地全文](src_dl/literature_2026-10-29_model/leechan2013/paper.pdf)。

- **核对**：Sec. II Eqs. (1)–(5)，Sec. III C Eqs. (28)–(31)，Sec. IV Eq. (32)、Fig. 1。逐原子独立衰减、平均场三次方程、线性涨落和相空间展开均已研究。
- **落点**：Model 平均场方程之前，及 Discussion 对局部线性涨落适用范围的说明。
- **边界**：线性化描述分支内涨落，不能自动解决分支间稀有切换或本文的激发模混合；作者也指出双稳态边界附近线性涨落展开失效。不能把平均场方程或通常的 Jacobian 分析列为本文新方法。

### 3. Roberts–Clerk 2023：同模型已有有限 N 精确稳态

**D. Roberts and A. A. Clerk, “Exact Solution of the Infinite-Range Dissipative Transverse-Field Ising Model,” Phys. Rev. Lett. 131, 190403 (2023).** [DOI](https://doi.org/10.1103/PhysRevLett.131.190403) · [arXiv](https://arxiv.org/abs/2307.06946) · [本地全文及补充材料](src_dl/literature_2026-10-29_model/roberts2023/paper.pdf)。

- **核对**：Eqs. (1)–(7)、(14)–(16)，Fig. 3，SM III.1。其模型族包含本文，精确稳态写作 ρ_ss ∝ ΦΦ†；不要求仅在最大自旋扇区，也能处理更一般的场分布。
- **落点**：Introduction 已知结果段；Model 模型来源；Discussion 稳态势和激发谱之间的区别。
- **创新性影响**：这是最重要的同模型遗漏之一。论文明确聚焦稳态，没有给出完整激发 Liouvillian 谱；因此不覆盖我们的谱避交叉，但必须承认稳态的解析基础已经存在。

### 4. Noachtar–Clerk 2026：同模型族的势垒—慢速率联系

**D. D. Noachtar and A. A. Clerk, “Exact metastability in a class of driven-dissipative quantum many-body systems,” arXiv:2606.07736v1 (2026).** [预印本](https://arxiv.org/abs/2606.07736) · [本地全文](src_dl/literature_2026-10-29_model/noachtar2026/paper.pdf)。截至本轮元数据核查，未列正式期刊出处。

- **核对**：Conjecture 1 / Eq. (37)，Sec. V B Eq. (48)、Eqs. (58)–(60)，Figs. 7–9。利用 hidden time-reversal symmetry 和特定稳态纯化定义势垒，提出并在例子中验证通常 Liouvillian gap 的指数标度。
- **必须说清**：动力学公式明确是 **conjecture**；对象是相对零模的最慢衰减率 Γ_diss = −max_{λ≠0} Re λ，而非本文 g_N = (λ₊−λ₋)/2。其自旋验证图含局域和集体衰减，不能直接写成已验证本文纯局域衰减的参数点。
- **落点**：Introduction 的 instanton/慢切换背景；Discussion 当前 Φ(n_f) 段。
- **影响**：不能笼统写“稳态信息不可能预测动力学”。应写：已有满足特定对称性、使用特殊纯化权重的预测；这些结果尚未建立本文激发数边缘率函数与有限率谱间距的关系。特殊纯化势也不等于 −log ρ_ss 或任意人口边缘势。

### 5. Xiang–Li–Bai–Ma 2026：最近的 Rydberg 切换与谱研究

**Y.-X. Xiang, W. Li, Z. Bai, and Y.-Q. Ma, “Switching metastable dynamics in many-body open quantum systems,” National Science Review 13, nwag146 (2026).** [正式全文](https://academic.oup.com/nsr/article/13/10/nwag146/8512294) · [DOI](https://doi.org/10.1093/nsr/nwag146) · [本地预印本](src_dl/literature_2026-10-29_model/xiang2025/paper.pdf)。预印本 arXiv:2505.05202 的题名不同；引用正式版题名。

- **核对**：正式版全文 Model Eqs. (1)–(5) 后的量子亚稳态小节与 Fig. 2；预印本 Secs. II C–D、Fig. 3。研究切换、初态依赖、通常 Liouvillian gap 及尺寸标度。
- **关键区别**：其集体计算采用 S=N/2 空间及 L_eff = Σ_M √(M+S)|M−1⟩⟨M|；正式文中明确称为 *heuristic jump operator*，忽略单原子跳跃造成的相位差。不能把它等同于本文保留所有 j 块及块间转移的 Σ_i D[σ_i⁻]。
- **落点**：Introduction 切换 gap 与激发模 splitting 的比较；Model 完整置换不变空间的范围说明。
- **出版信息**：Crossref 核对为卷 13、文章号 nwag146、2026 年；网页显示的月份/期号存在不直观组合，BibTeX 不采用月份以免传播歧义。

### 6. Rubio-García et al. 2022：有限非零衰减率处近简并已经出现过

**Á. Rubio-García et al., “Exceptional Spectral Phase in a Dissipative Collective Spin Model,” Phys. Rev. A 106, L010201 (2022).** [DOI](https://doi.org/10.1103/PhysRevA.106.L010201) · [arXiv](https://arxiv.org/abs/2202.09337) · [本地全文](src_dl/literature_2026-10-29_model/rubio2022/paper.pdf)。完整作者见 BibTeX。

- **核对**：Eqs. (2)–(4)，Fig. 2(a,b)，Eq. (10) 与 Fig. 3。模型 H=−hJ_z，集体 J₊、J₋、J_z 耗散。M=0 谱双重态在大尺寸下趋于 λ_{2n−1,0}=λ_{2n,0}=−2|p|Γn，公共衰减率有限且非零。
- **关键细节**：Fig. 2(b) 所展示的指数衰减量是归一化右本征矢距离 d₁；不能把这张图误引为已拟合本征值分裂的指数。Fig. 2(a) 展示谱对靠拢。
- **落点**：Introduction 的 Cabot 2024 之后；Discussion 有限率谱对的意义。
- **创新性影响**：排除“有限衰减率处近简并和 EP 前驱是新现象”的宽泛表述。本稿可区分的对象是独立局域耗散 Rydberg 模型中随失谐变化的避交叉、物理响应交换及其相干消元结构。

### 7. Molina 2026：相反大残数与初态可见性已有直接公式

**R. A. Molina, “Spectroscopic signatures of a Liouvillian exceptional spectral phase in a collective spin,” SciPost Phys. Core 9, 045 (2026).** [DOI](https://doi.org/10.21468/SciPostPhysCore.9.3.045) · [正式 PDF](https://scipost.org/10.21468/SciPostPhysCore.9.3.045/pdf) · [本地正式全文](src_dl/literature_2026-10-29_model/molina2026_published/paper.pdf)。2026-07-30 正式发表；不能继续仅列 arXiv:2602.01375。

- **核对**：Secs. 3–4；正式 PDF 第 18 页 Appendix A.8，Eqs. (A.52)–(A.57)。作者明确写出近简并模残数 R_± = ±B/(2δ)+A/2+O(δ)，其响应合成为 A/(z−λ₀)+B/(z−λ₀)²；时域对应 (a+bt)e^{λ₀t}。一般初态与稳态激发的可见性不同。
- **落点**：Physical response 的残数公式之后、相反大残数段；Dynamical reach 的制备依赖段。
- **创新性影响**：残数可很大、符号相反、相消后物理信号有界，以及 EP 信号依赖 source/detector，均不能作为本文独有机制。具体高分支乘积态扰动的响应交换和无拟合动力学误差仍是本文模型的定量结果。
- **适用边界**：其模型为集体自旋耗散；本文没有验证自身残数满足严格 1/g_N 渐近律，也没有实轴上的精确 EP。引用其公式作为背景，不把我们两组有限 N 数字升级为已证 Jordan 极限。该文对有限尺寸可对角化的某些宽泛措辞也不应照搬，其特例 p=1 另有精确 EP。

### 8. Li–Sun–Yi 2026：必须回应的易为方向后续

**C. Li, K. Sun, and W. Yi, “Critical exponents near the Liouvillian exceptional structure in Rydberg vapors,” Chin. Phys. B 35, 030302 (2026).** [DOI](https://doi.org/10.1088/1674-1056/adfb57) · [正式全文](https://cpb.iphy.ac.cn/en/article/pdf/preview/10.1088/1674-1056/adfb57.pdf) · [本地全文](src_dl/literature_2026-10-29_model/cpb2026_pdf/paper.pdf)。

- **核对**：Sec. 2 Eqs. (3)–(6)、Fig. 2；Sec. 3 Eqs. (10)–(23)、Fig. 3。其 1/2、1/3 指数是稳态激发密度根随参数扰动的标度；论文以 nonlinear exceptional structure 描述 stationary cubic 的根合并。
- **落点**：Introduction 解释三种对象的段落；Model 的 fold Jacobian 特征多项式之后。
- **表述建议**：先准确陈述其非线性稳态根结果，再说明本文另外计算平均场稳定性 Jacobian 和有限 N 线性多体生成元。它没有做本文的 Jacobian 缺陷检验，也没有计算有限 N 多体谱。不要仅因作者用了“EP”便将三者当成同一对象，也不宜不加限定地称其整篇结果错误。

## 支撑论证的理论与方法：优先新增 7 篇

| 编号、文献 | 已核对的直接依据 | 建议落点与边界 |
|---|---|---|
| **9. D. Finkelstein-Shapiro et al.**, *Adiabatic elimination and subspace evolution of open quantum systems*, PRA **101**, 042102 (2020). [DOI](https://doi.org/10.1103/PhysRevA.101.042102)；[全文](src_dl/literature_2026-10-29_methods/finkelstein2020/paper.pdf) | Eq. (7) 给出 L_eff(z)=PLP+PLQ(z−QLQ)⁻¹QLP；Eqs. (12)–(14) 给出谱参数导数对归一化的修正；Appendix F 给出 T(z)=zI−L_eff(z) 及 w†T′(λ)v=1。 | Introduction 方法定位、Sec. V A/C。比只引 Feshbach 1958 直接得多。其主应用是慢子空间及稳态归一化，不能替代本文有限率根的实际计算。 |
| **10. S. Güttel and F. Tisseur**, *The nonlinear eigenvalue problem*, Acta Numerica **26**, 1–94 (2017). [DOI](https://doi.org/10.1017/S0962492917000034)；[全文](src_dl/literature_2026-10-29_methods/guttel2017/paper.pdf) | Sec. 2.4，Theorem 2.8 / Eq. (2.12) 的 Keldysh 展开；简单根的归一化为 w*F′(λ)v=1；Sec. 2 的 Jordan chains、代数与几何重数。 | Sec. V C 的 WᵀR=lᵀ(I−K_z)r 及双根条件。我们的关系是这一标准导数配对在 Schur 回代中的具体表现，不能称新的一般定理。此处 Keldysh 是非线性谱理论，不是先前研究中的非平衡 Keldysh 路径积分。 |
| **11. J. L. García-Palacios and D. Zueco**, *Solving spin quantum master equations with matrix continued-fraction methods: application to superparamagnets*, J. Phys. A: Math. Gen. **39**, 13243–13284 (2006). [DOI](https://doi.org/10.1088/0305-4470/39/42/005)；[全文](src_dl/literature_2026-10-29_methods/garciapalacios2006/paper.pdf) | Sec. 7；Appendix G，Eqs. (G.1)–(G.5)：三项块递推、有限链从末端消元、Laplace 参数与矩阵连分式。 | Sec. V A 递推式与 SM coherence recursion。本文的 q=m−m′ 分组及保留全部 j 块是具体实现；矩阵连分式消元本身已有成熟方法。原文自旋-热浴模型与本文不同。 |
| **12. W. D. Heiss**, *The physics of exceptional points*, J. Phys. A: Math. Theor. **45**, 444016 (2012). [DOI](https://doi.org/10.1088/1751-8113/45/44/444016)；[全文](src_dl/literature_2026-10-29_methods/heiss2012/paper.pdf) | Sec. 2 Eqs. (1)–(12)：一般复参数 2×2 模型、耦合乘积、平方根分支、自正交、Jordan 结构及 resolvent 二阶极点。 | Sec. IV 两模判别式与复失谐 EP；Sec. V C。两圈复位只指本稿的本征值分支，不能据此宣称归一化本征矢相位也在两圈后复位。 |
| **13. V. V. Albert and L. Jiang**, *Symmetries and conserved quantities in Lindblad master equations*, PRA **89**, 022118 (2014). [DOI](https://doi.org/10.1103/PhysRevA.89.022118)；[全文](src_dl/literature_2026-10-29_methods/albert2014/paper.pdf) | Appendix B 在 Eq. (B3) 后说明零本征值没有广义本征矢；Appendix C 以有界演化排除虚轴上非平凡 Jordan 块。 | Introduction 的有限维 CPTP 半群零模不可缺陷化。这是既有结构定理，不是本文的新发现；不用于非物理复失谐延拓。 |
| **14. T. Mori and T. Shirai**, *Resolving a Discrepancy between Liouvillian Gap and Relaxation Time in Boundary-Dissipated Quantum Many-Body Systems*, PRL **125**, 230604 (2020). [DOI](https://doi.org/10.1103/PhysRevLett.125.230604)；[全文](src_dl/literature_2026-10-29_methods/mori2020/paper.pdf) | Eqs. (3)–(10)、(19)：展开系数和左右模小重叠可使有限时弛豫远慢于普通谱隙直觉。Eq. (10) 明写系数分母中的左右模配对。 | Sec. III 残数与 Sec. VI 直接动力学验证的动机。边界耗散晶格不是本文模型，不应照搬其系统尺寸指数。 |
| **15. T. Mori**, *Metastability associated with many-body explosion of eigenmode expansion coefficients*, PRResearch **3**, 043137 (2021). [DOI](https://doi.org/10.1103/PhysRevResearch.3.043137)；[全文](src_dl/literature_2026-10-29_methods/mori2021/paper.pdf) | Sec. II 的 L¹/L∞ 对偶归一化、Eqs. (4)、(8)–(10)，及双井相互作用粒子/FA 模型反例：大展开系数可以造成慢弛豫而不需要趋零本征值。 | Introduction 亚稳态定义的限定；Sec. VI/Discussion。研究对象是经典随机系统；用于反驳过宽的必要条件，不用于宣称本文已经发现同一种长寿命机制。 |

## 哪些常见表述需要随引用一起校准

1. **“亚稳态需要小本征值”不能无条件写。** 当前 Introduction 第二段的 “Metastability instead involves small but nonzero decay eigenvalues and a separation of relaxation scales” 过宽。建议改为：

   > In the spectral-separation framework, metastability is associated with low-lying decay modes separated from faster relaxation. More generally, large eigenmode-expansion coefficients can also prolong transients, so eigenvalue spacings alone do not determine relaxation times.

   前句引 Macieszczak 2016/2021，后句引 Mori–Shirai 2020 和 Mori 2021。本文自身的 g_N 小只说明两条有限率谱靠近，既不能单独证明，也不能单独否定任意制备下的长时间尺度。

2. **“精确消元”是计算准确性，不自动是方法创新。** K(z)、矩阵连分式、导数归一化分别已有直接先例。可主张的是它们在当前局域耗散 Rydberg 谱对中的应用和验证：同一 μ_N(z) 的两个根、冻结 K(z_c) 的失败、误差控制与物理响应。

3. **“有限衰减率处的谱双重态”已有先例。** 引 Rubio-García 2022 和 Molina 2026 后，创新性表述应强调本稿的具体模型、实参数避交叉及响应交换，不应依赖“激发谱此前无人研究”。

4. **大残数与抵消已有公式。** 当前 N=32、40 的有符号数值可保留，但本稿尚未证明残数的 1/g_N 渐近律，也未证明 N→∞ 实轴 EP。建议明确区分“与既有近 Jordan 响应机制相容”和“已经建立该极限”。

5. **Φ(n_f) 段需回应 Noachtar–Clerk。** 现有局限判断仍成立：本文没有建立 σ=Φ(n_f)。但“稳态边缘分布不够”不能被读成“任何稳态信息都不能预测任何速率”。特殊纯化势、普通人口边缘势，以及有限率激发模 splitting 是三个不同对象。

6. **Xie/Li/Sun/Yi 的非线性稳态根与本文 Jacobian 应明确分开。** Li–Sun–Yi 的 1/2、1/3 结果属于定常根分岔；本稿 fold 上 χ(λ)=λ[(λ+1)²+Ω²/2] 属于线性稳定性谱。用定义和计算区分，避免用同一个“EP”词掩盖问题。

## 有限尺寸模型映射：确认“同模型”不是靠题目相似

本稿 N_e=J_z+N/2。省略常数，

\[
H=\Omega J_x+\frac{V}{2N}J_z^2+
\left[-\Delta+\frac{V(N-1)}{2N}\right]J_z.
\]

四个可直接对应的文献参数如下；所有局域衰减率均取本稿 γ。

| 文献约定 | 与本稿严格匹配的参数 |
|---|---|
| Lee 2012：V_L Σ_{i<j} n_i n_j/(N−1) | V_L=V(N−1)/N；Δ_L=Δ；Ω_L=Ω |
| Lee–Chan 2013：V_L N_e²/(2N)−Δ_L N_e | V_L=V；Δ_L=Δ+V/(2N)；Ω_L=Ω。其 J 定义为 Σσ，阅读集体自旋公式时另有因子 2。 |
| Roberts–Clerk 2023：J_RC Σ_{i≠j} S_i^z S_j^z + h^z J_z + h^x J_x | J_RC=V/(2N)；h^z=−Δ+V(N−1)/(2N)；h^x=Ω；集体衰减 Γ_RC=0 |
| Noachtar–Clerk 2026 Eq. (48)：\bar J J_z²/N + Δ_NC J_z + Ω_NC J_x；集体衰减系数 \barΓ/N | \bar J=V/2；Δ_NC=−Δ+V(N−1)/(2N)；Ω_NC=Ω；\barΓ=0 |

O(1/N) 失谐项不可在有限尺寸谱比较中遗漏。Li–Sun–Yi 的**平均场方程**则由 Ω_CPB=Ω/2、Δ_CPB=−Δ、V_CPB=V、Γ=γ 对应；这不是声称其全文已计算同一有限 N 线性生成元。

## 条件新增：7 篇中按实际论述选用

| 文献 | 用途、核对程度和建议 |
|---|---|
| **Kessler et al.**, *Dissipative phase transition in a central spin system*, PRA **86**, 012116 (2012). [DOI](https://doi.org/10.1103/PhysRevA.86.012116)；[全文](src_dl/literature_2026-10-29_model/kessler2012/paper.pdf) | 已核对 Sec. IV B、Fig. 10 附近“avoided crossing due to mode mixing”。避免把有限尺寸 Liouvillian 避交叉当作首见；其临界处公共速率趋零，耗散中心自旋模型不同。建议引言补入，优先于再添综述。 |
| **Dutta–Zhang–Haque**, *Quantum Origin of Limit Cycles, Fixed Points, and Critical Slowing Down*, PRL **134**, 050407 (2025). [DOI](https://doi.org/10.1103/PhysRevLett.134.050407)；[全文及 SM](src_dl/literature_2026-10-29_model/dutta2025/paper.pdf) | 已核对 Fig. 2 及正文 Jacobian 谱梯级论述；SM S5.2、S7。适合解释分支内局部松弛谱，与另一吸引子的指数小切换模分开。不同的集体跳跃和吸引子结构，不能直接预测本文数值。 |
| **Macieszczak–Rose–Lesanovsky–Garrahan**, *Theory of classical metastability in open quantum systems*, PRResearch **3**, 033047 (2021). [DOI](https://doi.org/10.1103/PhysRevResearch.3.033047)；[全文](src_dl/literature_2026-10-29_methods/macieszczak2021/paper.pdf) | 已核对 Secs. II–III，尤其 Eqs. (3)–(8) 和 classicality 条件。支持何时投影系数能解释为亚稳相概率；本稿有符号模残数通常不具备这种解释。建议在相应限定展开时加入。 |
| **Mori–Shirai**, *Symmetrized Liouvillian Gap in Markovian Open Quantum Systems*, PRL **130**, 230404 (2023). [DOI](https://doi.org/10.1103/PhysRevLett.130.230404)；[全文](src_dl/literature_2026-10-29_methods/mori2023/paper.pdf) | 已核对 Eqs. (11)–(12)，给稳态自关联瞬态衰减上界。若不讨论严格上界，Mori 2020 已足够，不必重复。该文数值例子是相互作用量子点，不应误写为同一 Rydberg 模型。 |
| **Heiss–Sannino**, *Avoided level crossing and exceptional points*, J. Phys. A: Math. Gen. **23**, 1167–1178 (1990). [DOI](https://doi.org/10.1088/0305-4470/23/7/022)；[全文](src_dl/literature_2026-10-29_methods/heiss1990/paper.pdf) | 已核对 Secs. 1–2，提供避交叉与复参数分支点的历史来源。若希望同时引原始论文与一般非厄米公式，可与 Heiss 2012 并引；不据此声称任意高维非厄米避交叉均满足同一局部近似。 |
| **S. Nakajima**, *On Quantum Theory of Transport Phenomena: Steady Diffusion*, Prog. Theor. Phys. **20**, 948–959 (1958). [出版页](https://academic.oup.com/ptp/article/20/6/948/1930693)；DOI 10.1143/PTP.20.948 | 本轮核查出版社元数据与摘要，未逐式核对原文。当前保留时域记忆方程时，适合作为投影法历史出处。不能声称原文直接推导了本稿的 coherence-order 消元。 |
| **R. Zwanzig**, *Ensemble Method in the Theory of Irreversibility*, J. Chem. Phys. **33**, 1338–1341 (1960). [出版页](https://pubs.aip.org/jcp/article/33/5/1338/206419/Ensemble-Method-in-the-Theory-of-Irreversibility)；DOI 10.1063/1.1731409 | 本轮核查出版社元数据与摘要，正文受访问限制。与 Nakajima 并引于记忆核式。不要误配成 1961 年的 *Memory Effects in Irreversible Thermodynamics*。 |

另外已检查 Minganti et al., PRA **100**, 062131 (2019) 的摘要及相关正文，确立 Hamiltonian EP 与含量子跳跃的 Liouvillian EP 不应混用。当前稿件直接从完整 L 出发，若不再专门展开 no-jump 对比，不必为篇数新增它。Albert et al. PRX **6**, 041031 (2016) 的主要重心是渐近子空间与响应，现有直接文献已足够；Brody 的 biorthogonal quantum mechanics 也不适合作为 Liouvillian 有符号残数“概率解释”的依据。不同模型的 Ising 链、皮肤效应、一般非厄米综述不作为本轮优先新增。

## 现有 12 篇如何处理

| 已有引用 | 建议 |
|---|---|
| Carr 2013, PRL 111, 113901 | 保留实验背景；模型源头补 Lee 2012，不用实验论文承担理论模型出处。 |
| Casteels–Fazio–Ciuti 2017, PRA 95, 012128 | 保留一阶耗散相变和慢速率背景，明确其与本文有限率 splitting 不同。 |
| Minganti et al. 2018, PRA 98, 042118 | 保留一般 Liouvillian 相变谱理论。 |
| Xie et al. 2025, Science Bulletin 70, 3345–3350 | 保留研究动机，补 DOI [10.1016/j.scib.2025.08.051](https://doi.org/10.1016/j.scib.2025.08.051)。这是热 Rydberg 蒸气实验，不是本文理想全连接模型的精确有限 N 验证。 |
| Sun–Yi 2024, AAPPS Bulletin 34, 22 | 可保留方向综述，[DOI](https://doi.org/10.1007/s43673-024-00129-3)。不能代替 Li–Sun–Yi 2026 原始论文或 Heiss 的数学背景。 |
| Macieszczak et al. 2016, PRL 116, 240404 | 保留谱分离型亚稳态基础；适用条件需要明确。 |
| Ptaszyński–Chudak–Esposito 2026, arXiv:2604.15091 | 必须保留。已复查 v2 的 Fig. 3 和 Liouvillian gap 段，预测对象是最慢切换速率。它不推导本文有限率两模分裂。 |
| Chudak–Esposito–Ptaszyński 2026, arXiv:2610.21684 | 必须保留。已复查摘要、Sec. I–II 和方法边界：包括局域/集体衰减、非高斯项和非凸动量 Hamiltonian。与本稿失败的直接 gMAM 尝试相关，但不能当作已有 σ=0.14 预测。 |
| Krimer–Pletyukhov 2019, PRL 123, 110604 | 保留 few-mode/hysteresis 背景，[DOI](https://doi.org/10.1103/PhysRevLett.123.110604)。本轮没有重新逐式审计其整个几何约化。 |
| Cabot–Giorgi–Zambrini 2024, PRX Quantum 5, 030325 | 必须保留。[DOI](https://doi.org/10.1103/PRXQuantum.5.030325)。已复查全文的谱分解和 EP/phase-diffusion/activation 背景；是不同驱动振子中的 EP–分岔联系，不是本文的有限率实谱对。 |
| Feshbach 1958, Ann. Phys. 5, 357–390 | 保留历史来源，[DOI](https://doi.org/10.1016/0003-4916(58)90007-1)，并补直接 Liouvillian 和非线性本征问题文献。 |
| Shammah et al. 2018, PRA 98, 063815 | 必须保留实现来源；正文和 SM 中完整 Dicke 块与局域耗散的说明继续明确。 |

## 可直接实施的引用布局

以下定位按当前 main.tex；尚未写入。

| 位置 | 建议引用与论证任务 |
|---|---|
| Introduction 第一段（约第 18 行） | Lee 2012、Lee–Chan 2013、Carr 2013、Xie 2025：模型和实验动机。 |
| Introduction 第二段（约第 20 行） | Li–Sun–Yi 2026：非线性稳态根结构；Albert–Jiang 2014：有限维半群零模；Macieszczak + Mori：限定亚稳态叙述。 |
| Introduction 已知结果段（约第 22 行） | Roberts–Clerk 2023；Noachtar–Clerk 2026；Xiang 2026；保留两篇 instanton 文献。随后用 Rubio-García、Molina、Cabot 区分已知谱现象。避免把所有论文塞进一个无比较的长括号。 |
| Model（约第 35–59 行） | Lee 两篇、Roberts–Clerk；Shammah；需要时 Dutta 2025。 |
| Physical response（约第 90–112 行） | Mori–Shirai 2020；Molina 2026；若讨论概率解释，补 Macieszczak 2021。 |
| Two-mode/complex EP（约第 124–183 行） | Heiss 2012；可选 Heiss–Sannino 1990。 |
| Coherence elimination（约第 188–217 行） | García-Palacios–Zueco；Feshbach；Finkelstein-Shapiro；记忆核式加 Nakajima–Zwanzig。 |
| Derivative norm/EP condition（约第 260–276 行） | Güttel–Tisseur、Finkelstein-Shapiro Appendix F、Heiss。 |
| Dynamical reach（约第 280–282 行） | Molina、Mori–Shirai；引用后仍保留本文独立计算的两模/四模误差，不以文献替代验证。 |
| Discussion（约第 297 行） | Roberts–Clerk、Noachtar–Clerk、两篇 instanton：稳态势、切换指数、激发谱间距之间尚未建立的联系。 |

建议贡献句：

> We resolve a finite-size avoided crossing between two finite-decay modes in the locally damped, uniformly interacting Rydberg model, identify an exchange of high-branch response under physical product-state preparations, and recover the pair as two self-consistent roots of a coherence-resolved Schur operator. The analysis quantifies the dynamical reach of this spectral pair without identifying its finite-window splitting exponent with a switching action.

这句话不把通用投影公式、残数相消或 EP 的平方根结构当成新发现，也不声称已证明大 N 机制。最有价值的定位是：**把已有一般理论落实到一个实验动机明确、保留完整局域耗散的模型，给出可信的有限尺寸谱与物理响应结果。** 文献调研支持继续按 PRA 谱结构论文完善，不足以恢复原先“普适 fold–EP–instanton 指数”式 PRL 主张。

## 证据与交付范围

- “全文核对”指取得全文并读了与本稿主张直接相关的章节、方程或图，不表示逐页审稿。关键新颖性比较（Lee、Roberts、Rubio-García、Molina、Li–Sun–Yi、Noachtar、Xiang）已经进行原文检查。
- 方法核查使用下载全文；Nakajima、Zwanzig 明确只核对出版页和摘要。下载成功本身不计为实质核对。
- 2026 年预印本保留版本和预印本身份；Molina、Xiang 使用正式出版信息，不能仅沿用旧 arXiv 元数据。
- 本地资料保存在两个独立下载树 [模型文献](src_dl/literature_2026-10-29_model/) 与 [方法文献](src_dl/literature_2026-10-29_methods/)。报告和 BibTeX 可独立使用，第三方全文未作为报告提交内容。
- 未修改正文、补充材料、PDF、图、数值数据或原源码包；未运行新数值研究、未投稿、未推送。
