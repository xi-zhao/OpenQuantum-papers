# OpenQuantum-papers

量子物理科研稿件、计算代码、原始数据、图表、文献调研和验证记录的独立仓库。2026-10-29 从 OpenQuantum 产品仓库迁入八个研究目录；保留原有 `outputs/<项目>/` 路径，以便继续使用已有复现命令。

目录名中的 `prl` 是历史选题名称，不表示论文已经达到 PRL 标准、投稿或发表。下表反映迁移时已有稿件和记录的状态；迁移本身不构成新的科学验证。

## 稿件与研究记录

| 项目 | 当前入口 | 状态与边界 |
|---|---|---|
| Rydberg 系综中有限尺寸衰减模混合 | [PRA 稿件与复现说明](outputs/prl-chiral-switching/PRA_README.md) · [正文](outputs/prl-chiral-switching/paper/main.pdf) · [补充材料](outputs/prl-chiral-switching/paper/supplement.pdf) | PRA 草稿；已补充文献归属与有限尺寸限定，未证明间距指数等于作用量。 |
| 费米 Liouvillian EP 附近依可观测量而异的弛豫 | [研究说明](outputs/prl-highorder-ep/README.md) · [正文](outputs/prl-highorder-ep/paper/manuscript.pdf) · [补充材料](outputs/prl-highorder-ep/paper/supplement.pdf) | PRA 草稿；保留可观测量选择、条件化制备代价及半经典适用范围。 |
| 随机编码中经典与量子信息的共同恢复转变 | [研究说明](outputs/prl-distinguishability/README.md) · [正文](outputs/prl-distinguishability/paper.pdf) · [补充材料](outputs/prl-distinguishability/supplement.pdf) | 现有稿件与数值验证记录；旧两态稿和 Lean 证明单独保留，不冒充当前恢复定理的形式化证明。 |
| Rydberg 链中的 E7 谱双通道信号 | [正文](outputs/e7-rydberg-prl/paper/e7_rydberg.pdf) · [修订记录](outputs/e7-rydberg-prl/paper/REVISION_NOTES.md) | 有限链稿件；A4 收敛和更重粒子分辨仍有明确限制。 |
| 耗散 SYK 的慢模竞争 | [正文](outputs/prl-dissipative-syk/paper/main.pdf) · [LaTeX](outputs/prl-dissipative-syk/paper/main.tex) · [文献调研](outputs/prl-dissipative-syk/paper/literature_review.md) | 现有 PRE 草稿；文献报告的待办与科学结论未经本次迁移重新验收。 |
| 临界态首光子统计 | [研究记录](outputs/prl-critical-noise/README.md) · [研究笔记](outputs/prl-critical-noise/paper/research_note.pdf) | 创新性否决记录；核心结果与既有文献重合，不列作待投 PRL。 |

## 早期探索

- [绝热手性探索](outputs/prl-adiabatic-chirality/code/)：四个早期计算脚本，尚无完成稿件。
- [中性原子非厄米选题预检查](outputs/prl-nh-neutral/PRECHECK_A.md)：保留否决结论、计算和文献边界。

## 使用与历史

各项目自己的 README、修订记录及验证报告是复现入口。计算环境按项目管理；本仓库没有统一安装步骤，也不依赖 OpenQuantum 的宿主、插件运行时或模型配置。不要为了阅读论文运行全部探索脚本，部分脚本会覆盖已有数据或执行较大的计算。

原有七条研究分支经路径筛选后保存在 `history/openquantum/c/*`，共 27 个科研提交。提交作者、日期、说明及科研文件内容保留；去除产品文件后，提交 SHA 发生变化。原 SHA 与新 SHA 的对应关系、迁移清单和验证结果见 [迁移说明](MIGRATION.md) 与 [来源记录](provenance/)。新仓库的工作分支整合迁移时磁盘上的最新版本，包括此前尚未提交的研究文件。

构建缓存、下载的第三方文献和可重建 ZIP 完整保留在本地，但由 `.gitignore` 排除。正式稿件 PDF、图、代码、数据和研究记录纳入版本管理。未添加新的开放许可；第三方文献的权利不因存放于此而改变。

OpenQuantum 中原来的八个路径保留本地兼容符号链接；后续科研维护请直接在本仓库进行。教学资源库、公司监控表和服务产品集成的工程实验继续留在 OpenQuantum。

GitHub 远端：[xi-zhao/OpenQuantum-papers](https://github.com/xi-zhao/OpenQuantum-papers)，可见性为私有。默认分支为 `main`；七条研究历史分支保存在 `history/openquantum/c/*`。2026-10-29 完成私有远端备份；文献下载和缓存仍仅在本地，本次未投稿或公开发布。
