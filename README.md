# OpenQuantum-papers

`main` 仅保存论文 PDF（含补充材料）和本索引。LaTeX 源码、图件、计算代码、原始数据、验证记录与旧稿保存在研究分支。目录中的 `prl` 是历史名称，不表示投稿或发表状态。

**当前版本：2026-11-11。六篇正文及既有补充材料统一署名 OpenQuantum，PDF 作者属性同步。**

| 论文 | 正文 | 补充材料 | 去重外部文献数 |
|---|---|---|---:|
| Rydberg 系综中的有限尺寸衰减模混合 | [PDF · 9 页](outputs/prl-chiral-switching/paper/main.pdf) | [PDF · 13 页](outputs/prl-chiral-switching/paper/supplement.pdf) | 30 |
| 费米 Liouvillian EP 附近依可观测量而异的弛豫 | [PDF · 4 页](outputs/prl-highorder-ep/paper/manuscript.pdf) | [PDF · 9 页](outputs/prl-highorder-ep/paper/supplement.pdf) | 32 |
| 随机编码中经典与量子信息的共同恢复转变 | [PDF · 5 页](outputs/prl-distinguishability/paper.pdf) | [PDF · 8 页](outputs/prl-distinguishability/supplement.pdf) | 33 |
| Rydberg 链中 E7 谱的双通道信号 | [PDF · 4 页](outputs/e7-rydberg-prl/paper/e7_rydberg.pdf) | [PDF · 7 页](outputs/e7-rydberg-prl/paper/e7_rydberg_supplement.pdf) | 31 |
| 耗散 Majorana 系统中四时关联的能量选择 | [PDF · 4 页](outputs/prl-dissipative-syk/paper/main.pdf) | [PDF · 8 页](outputs/prl-dissipative-syk/paper/supplement.pdf) | 24 |
| QwenQPU：在语言模型标准层内执行量子电路模拟 | [PDF · 17 页](outputs/qwenqpu/paper/main.pdf) | 附录含于正文 PDF | 28（另含 2 条项目文档与数据，共 30 条） |

外部文献数不包含 Supplemental Material（SM）入口。随机码的 33 项为 32 篇论文和 Feller 的一本书，含 SM 共 34 条参考项；上一版索引将外部文献误计为 34，现已更正。

## 2026-11-11 更新

- 新增 QwenQPU 英文论文：17 页、5 幅图、5 张正文表和 2 张附录表，署名 OpenQuantum。
- 已核对 30 条参考文献和 47 处正文引用，修正 2 处来源对应问题并完善书目信息；Wilson 区间、McNemar 检验和保守配对差区间复算与记录一致。
- 本次加入的是参考文献核查后的 PDF；论文描述经典硬件上的有限容量量子电路模拟，不主张物理量子加速。arXiv 投稿尚未提交。
- LaTeX 源码、5 幅图、源码 ZIP 与引用核验记录保存在 [QwenQPU 研究分支](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/qwenqpu-paper-20261111/outputs/qwenqpu)。原有十份 PDF 保持原样；主分支现有十一份 PDF，共 88 页。

## 2026-11-09 更新

- Rydberg 补准固定失谐与尺寸范围，图注和轴名明确绝对相对误差；费米 EP 补准纯指数尾部的实单根与重叠条件，保留 Jordan 多项式例外；随机码逐项明确 SDP 验证的尺寸—系综组合。
- E7 澄清反射对称包络、开链质量比的分母和图例；Majorana 统一算符编号并定义能量分量记号。既有有限尺寸、资源成本、必要对照及失败边界均保留。
- 五稿正文均有 SM 入口，且正文参考表覆盖各自 SI 的全部外部文献。E7 将原仅列于 SI 的五篇纳入正文表，去重仍为 31 篇。E7 的 Vovrosh 等文献更新为 PRX Quantum 7, 043018 (2026)，Majorana 的 García-García 等文献更新为 Phys. Rev. D 110, 086010 (2024)。
- 十份 PDF 共 71 页，已完成编译、引用、字体嵌入、署名和版面核查。本次同步逐份核对 SHA256，使用已经验证的 PDF，没有重跑物理计算，不改变投稿或发表状态。

## 版本来源

上述五篇物理论文的十份 PDF 来源于本地研究分支 `c/paper-review-fixes-20261109`，来源提交为 `d873e073f4f2247512ebd710b512ebabc4be6c9d`。完整源稿、书目、修复前审查报告及本轮验证记录保存在该研究提交中；2026-11-09 的同步仅将十份 PDF 和该版索引放入远端 `main`。

E7 的科学基线为研究分支 `c/e7-manuscript-evidence` 的 `085fc08c79917ba007fd08650b1acf3ee768c4a9`；后续修订保留约 62% 的绑定偏差、固定能隙检验未过门槛和重粒子未分辨等结论边界。其余稿件采用各自已验证版本。

研究分支中的 `outputs/paper-review-fixes-2026-11-09/README.md`、`validation.json` 和 `delivery-validation.json` 保存修复清单、引用口径、验证范围、实际构建提示和逐文件 SHA256。旧审计及旧包保留历史语义。本地排版审阅包不是完整数值复现包；原始样本和计算材料仍保存在研究资料中。

## 历史研究资料

- [完整迁移快照](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/import-research-papers)：保留迁移时的来源记录、研究笔记、失败尝试、代码和数据。
- [Rydberg 衰减模混合历史研究分支](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/chiral-pra-concise)：保留此前稿件及配套计算、验证和交付清单。

研究笔记和创新性否决记录不列为论文成稿。更新 `main` 时只同步需要交付的 PDF，不合并整个研究分支。下载的第三方文献与编译缓存不进入 `main`。
