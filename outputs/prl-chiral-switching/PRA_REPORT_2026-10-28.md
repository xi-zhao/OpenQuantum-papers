# PRA 稿件重写与验证记录

日期：2026-10-28。目标期刊：Physical Review A。

题目：**Finite-size decay-mode hybridization near a mean-field fold in an open Rydberg ensemble**

## 本轮交付与论证主线

正文和补充材料已按有限尺寸谱结构重新组织。主问题是：平均场折点附近的两个有限衰减实模如何交换物理响应，非厄米耦合如何控制它们的分裂，以及相干消元在什么条件下保留这对根。

论文只使用 Ω=2γ、V=10γ、独立局域衰减的统一双精度证据；完整生成元指完整置换不变算符空间，包含所有 Dicke j 块。没有把旧的多驱动经验拟合、未收敛瞬子、手性标度或预注册通过次数搬入当前主结论。

保留的结果：

- N=24,28,…,56 的实模 avoided crossing，有限窗口半间隙斜率 σ=0.14050；明确未确定 N→∞ 极限。
- 合法乘积态扰动的谱权重在交叉两侧交换高支响应；N=32 扫描与 N=40 对照一致。
- 双正交投影给出规范不变的 g²=V_ab V_ba；单向耦合指数依赖归一化，不能当成唯一物理作用量。
- 相干阶数消元给出依赖谱参数的 K(z)，两条衰减谱由同一有效支的两个自洽根恢复；冻结 K(z_c) 会丢掉这对邻近根。
- N=24、32、40 的复失谐 EP2 局部验证，及 N=32 一圈交换、两圈复位。
- 两实模不足以描述所测瞬态；加入更慢的振荡共轭对显著降低误差。

**没有建立的结论**：σ=Φ(n_f)、独立大 N 作用量预测、Δ_*(N) 必然趋于折点、普适 ghost-escape 模身份、短程相互作用推广、普适手性标度律。正文对此均有限定。

## 新增验证

| 验证 | 结果 | 证据 |
|---|---|---|
| 平均场 fold Jacobian | 在 F=F_n=0 下符号化简为 λ[(λ+1)²+Ω²/2]；零根简单 | `data/pra_spectral/meanfield.json` |
| 初态的概率与相干约定 | N=1,2,4 二项人口和单原子相干误差约 10⁻¹⁶ | `modes_N32.json` |
| 独立全自旋演化 | N=1–4、四个时刻，激发密度误差 <8.7×10⁻¹⁶，纯度误差 <1.5×10⁻¹⁵ | `validation_repeat.json` |
| 高支谱权重交换 | N=32 从 λ+ 主导变为 λ− 主导；N=40 相同 | `modes_N32.json`、`modes_N40_quick.json` |
| 初态步长减半 | 归一化谱权重变化 <7.3×10⁻⁶ | `validation_repeat.json` |
| 种子及谱位移改变 | 三点目标谱差 <9×10⁻¹⁵ | `validation_repeat.json` |
| 物理瞬态，无振幅拟合 | N=32 交叉处、t=6–12，高支误差 25.1%→0.506%；折点初态 69.5%→5.00% | `modes_N32.json` |
| 复 EP 与分支绕行 | 三尺寸及加密/换种子控制均通过；局部证据，不是全谱枚举 | `ep_summary.json` |
| 相干消元八个根 | 完整生成元回代最大残差 1.133×10⁻¹³；半间隙与独立谱相对差 <1.3×10⁻⁹ | `data/coherence_mechanism/validation.json` |
| 二次近似 | 半间隙相对误差从 N=24 的 1.20% 降到 N=48 的 6.59×10⁻⁶ | `fraction_derivatives_N*.json` |

本轮新增 `pra_validate.py` 将原先 inline 的步长/种子检查保存为可复现入口，同时用显式张量积全自旋构造和 DOP853 积分检查物理初态、Dicke 重数以及动力学。初版参考使用 dense `expm` 时触发本机 Accelerate 的浮点 warning；最终改用独立 ODE 参考，`python3 -W error ... --repeat-modes` 通过，未屏蔽 warning。

主会话核对了投影规范变换、Q 空间二阶修正、相干递推导数、WᵀR=−f_z 的根上条件，以及双根对应 EP2 所需的正则性和简单有效支条件。图注已明确：绝对人口剖面比例虽与模归一化无关，但不是物理态的概率；小谱间隙不表示该实模对有长寿命。

## 文献、审阅与限制

已复用最接近工作的下载全文，保留 Ptaszyński–Chudak–Esposito、Chudak–Esposito–Ptaszyński、Cabot–Giorgi–Zambrini 等引用，并明确 Schur/Feshbach 消元属于已有方法。补齐 Xie、Sun–Yi、Krimer–Pletyukhov 的完整文献题名。来源链接见 `PRA_README.md`。

最终独立审阅尝试因模型服务返回月度配额耗尽错误而提前终止，未产生审阅结论，也未擅自更换模型或 Provider。这里交付的是经过本地数学、数值、数据与排版核查的完整作者审读稿；不能写成“已通过独立审稿”，也不能保证 PRA 接收。

作者与机构仍保留显式占位。没有编造作者、基金、伦理或数据公开声明。源码包随稿本地交付，不等于已上传公开仓库或投稿。

## 工件与复现

最终正文 8 页、补充材料 6 页，包含五幅正文图和一幅补充图。两个 PDF 均经文字提取和逐页缩略图检查；所有字体嵌入，没有未定义引用或公式越栏。保留的 `nameref`/REVTeX 标签警告及 `h`→`ht` 浮动建议不影响显示。

- `paper/main.tex` / `main.pdf`：正文。
- `paper/supplement.tex` / `supplement.pdf`：推导和验证细节。
- `figs/pra_*.pdf` / `.png`：五幅正文图和一幅补充图。
- `PRA_README.md`、`requirements-pra.txt`：命令、版本、来源和证据索引。
- `data/pra_spectral/artifact_manifest.json`：包内文件 SHA256 与证据对应关系。
- `deliverables/pra-decay-mode-source-2026-10-28.zip`：明确文件白名单生成的源码与数据包，含 PDF、LaTeX、图、代码、JSON 和数值数组。

旧 PRL 正文/补充材料的 tex 与 PDF 均保留为 `paper/*_prl_before_pra_2026-10-28.*`，不混入当前源码包。失败瞬子及历史研究记录原位保留。

源码包已在 Elements 上的新临时目录解包，逐项核对 SHA256，然后运行数值一致性及小 N 对照、重新生成六幅图、编译两个 PDF，均通过；没有使用原研究目录中未打包的文件。生成的 ZIP 可由脚本重建，不重复纳入 Git，代码、数据和 PDF 纳入独立提交。

新稿在任务分支 `c/pra-decay-mode-hybridization` 交付，独立本地提交仅包含本轮文件。未 push、未创建 PR、未部署、未投稿。准确提交号以最终交付消息及 Git 记录为准。
