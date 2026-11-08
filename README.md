# OpenQuantum-papers

`main` 仅保存论文 PDF（含补充材料）和本索引。LaTeX 源码、图件、计算代码、原始数据、验证记录与旧稿保存在研究分支。目录中的 `prl` 是历史名称，不表示投稿或发表状态。

**当前版本：2026-11-08。五篇正文及全部补充材料统一署名 OpenQuantum，PDF 作者属性同步。**

| 论文 | 正文 | 补充材料 | 去重外部文献数 |
|---|---|---|---:|
| Rydberg 系综中的有限尺寸衰减模混合 | [PDF · 9 页](outputs/prl-chiral-switching/paper/main.pdf) | [PDF · 13 页](outputs/prl-chiral-switching/paper/supplement.pdf) | 30 |
| 费米 Liouvillian EP 附近依可观测量而异的弛豫 | [PDF · 4 页](outputs/prl-highorder-ep/paper/manuscript.pdf) | [PDF · 9 页](outputs/prl-highorder-ep/paper/supplement.pdf) | 32 |
| 随机编码中经典与量子信息的共同恢复转变 | [PDF · 5 页](outputs/prl-distinguishability/paper.pdf) | [PDF · 8 页](outputs/prl-distinguishability/supplement.pdf) | 34 |
| Rydberg 链中 E7 谱的双通道信号 | [PDF · 4 页](outputs/e7-rydberg-prl/paper/e7_rydberg.pdf) | [PDF · 7 页](outputs/e7-rydberg-prl/paper/e7_rydberg_supplement.pdf) | 31 |
| 耗散 Majorana 系统中四时关联的能量选择 | [PDF · 4 页](outputs/prl-dissipative-syk/paper/main.pdf) | [PDF · 8 页](outputs/prl-dissipative-syk/paper/supplement.pdf) | 24 |

## 本次更新

- 五篇稿件围绕各自最强的证据组织贡献叙事，保留有限尺寸、适用条件、必要对照与失败边界。
- E7 采用较新的正文＋补充材料结构，正文引用 26 篇、SI 引用 21 篇，去重后 31 篇；固定能隙检验及其有限尺寸结论完整保留。
- Majorana 正文新增 10 篇外部文献，目前为 24 篇外部文献及 1 条 Supplemental Material 入口；SI 引用 18 篇。新增引用分别对应费米回归、Markovian 生成元、耗散算符尺寸、有限时间相关界和多时相关测量，不以引文数量替代新颖性判断。
- 十份 PDF 均重新编译并核查作者属性、引用解析、字体嵌入及文件一致性；本次没有重跑物理计算，不改变投稿或发表状态。

## 版本来源

全部 PDF 来源于本地研究分支 `c/paper-polish-20261108`，来源提交为 `2acfd917263b3b031183090b2d8ec41e34782b47`。完整源稿、书目及本轮核查记录保存在该研究提交中；本次仅同步 PDF 和索引到远端 `main`。

E7 的科学基线为研究分支 `c/e7-manuscript-evidence` 的 `085fc08c79917ba007fd08650b1acf3ee768c4a9`，整合时保留其 31 篇文献及固定能隙检验，不以较早的独篇稿覆盖更新证据。其余稿件采用各自当前已验证版本。

研究分支中的 `outputs/paper-polish-2026-11-08/MAIN_SYNC_2026-11-08.md` 与 `MAIN_SYNC_VALIDATION.json` 保存引用计数口径、阅读范围、构建提示和逐文件 SHA256。下载的第三方文献与编译缓存不进入 `main`。

## 历史研究资料

- [完整迁移快照](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/import-research-papers)：保留迁移时的来源记录、研究笔记、失败尝试、代码和数据。
- [Rydberg 衰减模混合历史研究分支](https://github.com/xi-zhao/OpenQuantum-papers/tree/c/chiral-pra-concise)：保留此前稿件及配套计算、验证和交付清单。

研究笔记和创新性否决记录不列为论文成稿。更新 `main` 时只同步需要交付的 PDF，不合并整个研究分支。
