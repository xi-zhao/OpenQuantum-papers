# OpenQuantum-papers

`main` 仅保存论文 PDF（含补充材料）和本索引。LaTeX 源码、图件、计算代码、原始数据、验证记录与旧稿保存在研究分支。目录中的 `prl` 是历史名称，不表示投稿或发表状态。

| 论文 | 正文 | 补充材料 |
|---|---|---|
| Rydberg 系综中的有限尺寸衰减模混合 | [PDF](outputs/prl-chiral-switching/paper/main.pdf) | [PDF](outputs/prl-chiral-switching/paper/supplement.pdf) |
| 费米 Liouvillian EP 附近依可观测量而异的弛豫 | [PDF](outputs/prl-highorder-ep/paper/manuscript.pdf) | [PDF](outputs/prl-highorder-ep/paper/supplement.pdf) |
| 随机编码中经典与量子信息的共同恢复转变 | [PDF](outputs/prl-distinguishability/paper.pdf) | [PDF](outputs/prl-distinguishability/supplement.pdf) |
| Rydberg 链中 E7 谱的双通道信号 | [PDF](outputs/e7-rydberg-prl/paper/e7_rydberg.pdf) | [PDF](outputs/e7-rydberg-prl/paper/e7_rydberg_supplement.pdf) |
| 耗散 Majorana 系统中四时关联的能量选择 | [PDF](outputs/prl-dissipative-syk/paper/main.pdf) | [PDF](outputs/prl-dissipative-syk/paper/supplement.pdf) |

## 版本来源

- Rydberg 衰减模混合稿采用 [6cbf712](https://github.com/xi-zhao/OpenQuantum-papers/commit/6cbf712162f46bf4ec66b4526c9f870aaa3fe826)：署名 OpenQuantum，首页不显示日期。
- E7 双通道稿同步自研究分支 `c/e7-manuscript-evidence` 的 `085fc08`：正文四页、补充材料七页，署名 OpenQuantum，已删除单位编号、单位、日期和致谢占位内容。
- 耗散 Majorana 能量选择稿同步自本地研究分支 `c/syk-pra-energy-selection` 的 `36b3950`：PRA 定位，物理正文四页、数学补充材料七页，作者和单位仍为占位符；复现与验证记录保存在该研究提交中。
- 其余 PDF 沿用整理前 `main` 的 [d74b2f5 快照](https://github.com/xi-zhao/OpenQuantum-papers/commit/d74b2f5d67ed9c49be392e3224099a876ff0e00d)。本次整理不修改这些稿件，也不改变其验证或投稿状态。

## 研究资料

- [完整研究快照](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/import-research-papers)：保留全部项目的来源记录、研究笔记、失败尝试、代码和数据。
- [Rydberg 衰减模混合研究分支](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/chiral-pra-concise)：包含本篇最新稿件及配套计算、验证和交付清单。
- 耗散 Majorana 能量选择的源码、图表、计算代码与验证记录保存在本地研究分支 `c/syk-pra-energy-selection`，来源提交为 `36b3950080dbbe3068d7f3ad88959104deb4b229`。

研究笔记和创新性否决记录不列为论文成稿，可从研究分支查看。更新 `main` 时仅同步需要交付的 PDF，不合并整个研究分支。
