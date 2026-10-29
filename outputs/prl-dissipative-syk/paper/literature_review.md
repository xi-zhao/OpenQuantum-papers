# 文献调研：缺引文献与创新性评估

日期：2026-10-27。对象：`paper/main.tex`（PRE 版，《Energy mode versus two-fermion continuum》）。

## 调研范围与方法

- **arXiv API 全量检索**：截至 2026-10，检索 SYK 与 Lindblad / Liouvillian / Markovian 的组合，得到 25 篇开放 SYK 论文。
- **全文核查**：下载其中约 50 篇 PDF，用 `pdftotext` 转成文本后，grep 关键公式和关键词，包括 4μ、qμ、Nμ/2、ladder、Bethe、bound state、spectral gap、conservation law 等。
- **三条检索线并行**：
  1. Lindbladian SYK 的谱隙文献；
  2. 开放或双站点 SYK 的四点函数与 ladder；
  3. 一般 Liouvillian 谱隙理论，包括谱界、守恒量和能量流。
- **抽查**：对最关键的几条结论，我直接在全文中复核了原句或原式，下文标注“已核全文”。

未覆盖：
- 监测 SYK 的测量诱导相变（MIPT）系列，例如 Jian–Liu–Chen–Swingle–Zhang，它们用复本 ladder。
- Brownian SYK 的 Lindbladian 文献。
- 混沌系统能量扩散的经典文献，包括 Wilkinson、Ott、Jarzynski。

这些方向只做了检索，没有读全文。

---

## 一、必须补引（直接相关，审稿人会查）

| 文献 | 与本文的关系 | 核查 |
|---|---|---|
| Yoshimura & Sá, *Nat. Commun.* **15**, 9808 (2024), arXiv:2312.00649 | 原文：“the presence of conservation laws would spoil anomalous relaxation… We will report a systematic study of the role of conservation laws in anomalous relaxation in forthcoming work.” 他们定性提出、并留作后续的问题，正是本文定量回答的问题。**这是最重要的定位依据。** | 已核全文 |
| Barthel & Zhang, *PRA* **105**, 052224 (2022), arXiv:2012.05505 | Eq. (21)：μ_max ≥ Re λ_k ≥ μ_min，界由 Liouvillian 的 Hermitian 部分给出。本文的界 2μ ≤ Δ 是它的直接特例，应写成“应用已知结果”，不能作为新定理。 | 已核全文 |
| Cai & Barthel, *PRL* **111**, 150403 (2013), arXiv:1304.6890 | 开放 XXZ 链中，守恒的磁化使退相干时间在热力学极限下发散。这是本文“更低权重的守恒量胜过能量模”这一判断的标准前例。 | 已核全文 |
| Mori & Shirai, *PRL* **125**, 230604 (2020), arXiv:2006.10953 | 谱隙与实际弛豫时间的区别，属于标准参考。注意卷号是 125，不是 124。 | 已核摘要页 |
| Ozaki, Yoshida, Katsura, arXiv:2606.08079 (2026) | 与本文同一模型（线性跳跃），但哈密顿量不含无序。精确求解后发现，两点函数的谱隙在先取 N→∞ 时有限，可积系统也是如此。他们的谱隙只取对两点函数有贡献的本征值，见其 Eq. (73)。 | 已核全文 |
| Kawabata, Kulkarni, Li, Numasawa, Ryu, *PRB* **108**, 075110 (2023), arXiv:2210.04093 | 同一模型。研究耗散形式因子中的动力学相变，其中 μ≈0.35–0.37 的转变需要和本文的 μ_c 区分开。 | 已核全文 |
| Nosaka & Numasawa, arXiv:2009.10759 | 双耦合 SYK 的 OTOC ladder，按对称和反对称扇区分解，结构上与本文的交换扇区一致，应在方法部分引用。 | 已核全文 |
| Pelliconi, Lapierre, Ryu, arXiv:2608.19310 (2026) | Lindblad Yukawa-SYK 的 OTOC Bethe–Salpeter 方程，用束缚态定 λ_L。属于技术前例。 | 已核全文 |
| Liu & García-García, arXiv:2608.27897 (2026) | 同一课题组，在 Keldysh 轮廓上推导一般四点函数的 ladder 递推（其 SM Eq. S74–S78）。但只投影到 OTOC 轮廓上做计算。这是本文 ladder 方法最接近的前作，必须引用。 | 已核全文 |

## 二、建议补引（背景，增强可信度）

- **Macieszczak, Guţă, Lesanovsky, Garrahan**, *PRL* **116**, 240404 (2016)：亚稳态与慢流形投影，支撑本文的慢扇区/快扇区划分。
- **Žnidarič**, *PRE* **92**, 042143 (2015)：耗散多体系统的弛豫时间。
- **Mori & Shirai**, *PRL* **130**, 230404 (2023)：对称化 Liouvillian 谱隙。
- **Zhang & Yu**, *PRL* **130**, 250401 (2023)：全连接模型加热库时，耗散与 scrambling 之间的算符尺寸转变。
- **Denisov et al.**, *PRL* **123**, 140403 (2019)；**Sá, Ribeiro, Prosen**, *PRX* **10**, 021019 (2020)；**Can et al.**, *PRL* **123**, 234103 (2019)：随机 Liouvillian 谱的背景文献。
- **Almheiri, Milekhin, Swingle**, arXiv:1912.04912：SYK 与热库之间的能量流约束，和本文能量模的物理图像相关。
- **Maldacena & Milekhin**, *JHEP* **04** (2021) 258：实时虫洞形成；**Chen, Zhai, Zhang**, *JHEP* **07** (2017) 150：SYK 耦合热库。
- **Temme**, *J. Math. Phys.* **54**, 122110 (2013)：Davies 生成元的谱隙下界。
- **Kawabata et al.**, *PRX Quantum* **4**, 030328：耗散混沌的对称分类。可选。

---

## 三、需要正面回应的三处部分重叠

### 1. García-García, Sá, Verbaarschot, Zheng, *PRD* **107**, 106006 (2023)，已引用

- Eq. (E8) 写的是 “for small µ and sufficiently large N (N > 8)”，第一激发态为 E₀+4μ。所以 **N≤8 的例外他们已经注意到**，只是没有解释。本文的贡献是给出原因：N=6 时 PH 模取到 2μ；N=8 时 Nμ/2 恰好等于 qμ。文中应写成“解释了 [Garcia2023] 注意到的 N>8 条件”。
- 他们明确写了 “progress… requires… a more precise understanding of the relation between the spectral gap Δ and the decay rate E_g”。**本文正是在回答这个公开问题**，摘要和引言应直接引这句话。
- 他们的 μ_c≈0.15 是 Green 函数振荡消失的位置，即例外点，不是谱隙交叉，需要在文中区分。本文已写到 μ_EP≈0.16 并加以区分，可再引原文强化。

### 2. Zheng, Dukelsky, Molina, García-García, *PRD* **113**, 066018 (2026)，已引用

- 他们的谱隙 Γ₀ 是在 “the block with opposite parity to the steady state block” 中计算的，也就是**奇宇称 block（费米子相干）**，不是物理密度矩阵所在的偶扇区。
- 他们发现 “intruder” 实本征值与 SYK 本征态有关，并认为这些模在热力学极限下 “lose importance”。**这和本文的结论表面上矛盾**，原因是两边讨论的不是同一个扇区：他们只看两点函数所在的扇区，那里的 intruder 模确实不重要；本文讨论物理偶扇区，H 是其中的精确模，并且决定谱隙。
- 他们在 μ~0.3 看到 Loschmidt 回波的一阶转变变成 crossover，数值上接近本文的 0.323，但他们归因于例外点大量出现。文中要说明这两件事是否有关联；目前的证据不足以判断，应如实写“数值接近，机制未定”。

### 3. García-García, Lu, Sá, Verbaarschot, *PRE* **113**, 024116 (2026)，已引用

- H、P、PH 三个精确模出自这篇。全文中 grep 不到任何把 H 与谱隙联系起来的表述，所以“H 何时决定谱隙”这个问题不是他们提出的。
- 本文 N=6 的例外正来自他们的 PH 模，应写明出处。

---

## 四、创新性评估

### 已有前作，不能声称是本文的贡献

| 内容 | 出处 |
|---|---|
| L(H) = −qμH 精确，以及 P、PH 模 | PRE 113, 024116 (2026) |
| 有限 N 下 4μ 能隙，且 N>8 | PRD 107, 106006 (2023)，Eqs. (E8)、(E27) |
| μ→0 时衰减率有限（反常弛豫）及其 Keldysh 虫洞解释 | PRD 107, 106006 (2023) |
| 守恒律破坏反常弛豫（定性） | Yoshimura & Sá (2024) |
| 守恒量是耗散子本征算符，因而严格按指数衰减 | Teretenkov & Lychkovskiy (2024) |
| Re λ 的界由 Hermitian 部分给出 | Barthel & Zhang (2022) |
| 低权重守恒量使谱隙闭合 | Cai & Barthel (2013) |
| 开放 SYK 的 ladder 方法（OTOC 扇区） | 2403.12359、2608.27897、2608.19310 |
| 大 q 的 Γ₁ 公式 | Kulkarni et al. (2022) Eq. (22) |

### 目前检索范围内没有找到前作，可作为本文贡献

1. **时序（非 OTOC）扇区的四点函数 ladder，并用它求偶扇区的 Liouvillian 谱。** 三条检索线都确认，前作的 ladder 只用于 OTOC。
2. **在大 N ladder 中复现能量模束缚态 E=qμ。** 这是精确模与 Keldysh 虫洞（大 N）图像之间的第一个直接联系，也独立检验了 kernel 和符号约定。
3. **数值排除低于 2Γ₁ 的其他集体模**，单态、非单态和复数 E 都检查过。这是数值结论，不是证明。
4. **谱隙判据 qμ_c = 2Γ₁(μ_c)**，q=4 时 μ_c=0.323J，q=8 时 μ_c=0.0199J，大 q 下有闭式 4𝒥/(q√(q²−4))；以及 μ_c 以上谱隙是连续谱边缘、不是孤立本征值这一刻画。
5. **慢扇区/快扇区的系统划分**，包括：相干项一阶速率精确为 Nμ/2；N≥10 时 H 是最慢的 secular 模；对 N≤8 例外的解释；带模式识别的有限 N 交叉 μ_c(N)。

### 判断

- 创新点集中在**方法**（时序 ladder）和**定量判据**（μ_c 及其机制）。核心物理图像，也就是“守恒量会改变反常弛豫”，已有定性前作。
- 本文最好的定位是：**回答 PRD 107, 106006 和 Yoshimura–Sá 两篇明确留下的开放问题**。这个定位有据可查，审稿人也容易认可。
- 定量上的薄弱点：有限 N 的 μ_c 外推只能夹住 0.323，不能确认它；“没有其他束缚态”是数值结论。在 PRE 的标准下这两点可以接受，但都必须写明。
- 发表级别：PRE 或 PRB 合适；仅凭现有结果，不足以投 PRL 或 PRX。

---

## 五、建议对稿件的具体修改（未执行）

1. **摘要和引言**：引用 PRD 107 关于 Δ 与 E_g 关系的原句，以及 Yoshimura–Sá 关于守恒律的原句，把问题表述为回答这两个公开问题。
2. **第 II 节**：界 2μ ≤ Δ 改写为 Barthel–Zhang Eq. (21) 的应用；N≤8 的段落写明这是在解释 PRD 107 的 “N>8” 条件；XXZ 的例子补引 Cai–Barthel。
3. **第 IV 节**：ladder 方法补引 2608.27897（Keldysh 轮廓上的一般四点递推）、2009.10759（交换扇区分解）和 2608.19310。
4. **新增一段**，与 Zheng et al. 2510.15793 比较：说明两边的扇区不同，intruder 模的作用不矛盾，并讨论 μ~0.3 的数值接近。
5. **讨论部分**：补引 Ozaki–Yoshida–Katsura（可积系统也有反常弛豫）和 Kawabata et al. 2210.04093（耗散形式因子相变），并与本文的 μ_c 区分。
6. **参考文献**：按第一、二节补全，并核对卷号。
