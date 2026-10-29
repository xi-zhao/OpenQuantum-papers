# OpenQuantum 科研成果迁移记录

日期：2026-10-29。来源仓库为 `xi-zhao/OpenQuantum`；迁移起点为分支 `c/pra-decay-mode-hybridization`，提交 `6dd02de9a9fbce4a56c7833170b4dabad6d88a1b`。迁入本地同级目录 `OpenQuantum-papers`，不配置远端或发布内容。

## 范围与工作树

完整移动八个目录，共 1,578 个普通文件、226,527,794 字节（约 216 MiB）。保留目录名与内部层级，未改写任何原有科研文件。迁移时 `prl-highorder-ep` 的 11 个已跟踪未提交修改，以及各项目尚未跟踪的文件，均以磁盘上的实际版本保留。

| 目录 | 普通文件数 | 字节数 |
|---|---:|---:|
| `outputs/e7-rydberg-prl` | 179 | 19,463,077 |
| `outputs/prl-chiral-switching` | 723 | 157,236,470 |
| `outputs/prl-critical-noise` | 72 | 15,484,235 |
| `outputs/prl-dissipative-syk` | 246 | 2,227,108 |
| `outputs/prl-distinguishability` | 251 | 26,466,435 |
| `outputs/prl-highorder-ep` | 94 | 5,602,087 |
| `outputs/prl-adiabatic-chirality` | 4 | 4,153 |
| `outputs/prl-nh-neutral` | 9 | 44,229 |

其中 758 个科研文件纳入当前 Git 版本；其余 820 个文献下载、缓存和生成物继续保留本地，由忽略规则排除。原本已经版本化的历史稿 `archive/initial-draft/paper.bbl` 作为例外保留，避免因通用构建忽略规则遗漏旧版书目。否决选题、失败尝试和旧稿一并保留，状态见根 README 与各项目记录。

采用同一文件系统内的目录重命名，文件 inode、模式、修改时间、硬链接计数和扩展属性均保留。迁移前后逐文件 SHA256 一致。已有的 `prl-distinguishability/lean/.lake` 是指向 Elements 上 Lean 缓存的符号链接，链接本身完整保留且仍可解析；缓存不是本仓库内容，也不纳入 Git。

原 OpenQuantum 中的八个路径留下本地兼容符号链接；原仓库的 `.gitignore` 排除这些链接，其任务分支提交删除旧的跟踪文件并添加迁移说明。兼容链接仅服务本机旧路径，不随产品仓库克隆分发。后续科研工作应在本仓库进行，避免在原仓库切回旧科研分支而恢复第二份副本。

教学资源库、公司监控表、产品集成实验及产品代码保持在 OpenQuantum；未迁入产品凭据、依赖环境或运行时配置。

## 历史保留

使用 Git 自带的 `fast-export` / `fast-import`，只筛选上述八个路径，不安装额外过滤工具、不改写原仓库。导入全部七条相关本地分支和 27 个科研提交，分支保存在 `history/openquantum/c/*`。

| 原分支后缀（`c/…`） | 原分支头 | 导入后的分支头 |
|---|---|---|
| `chiral-biorthogonal-sigma` | `7f0f8be1e0e3` | `c66b79c9ebae` |
| `distinguishability-rigorous-revision` | `18615df8cad0` | `8e86fdfc819d` |
| `pra-decay-mode-hybridization` | `6dd02de9a9fb` | `542f9b790882` |
| `pra-ep-scaling` | `58655de1c2f0` | `f7ebccf406e1` |
| `pra-observable-relaxation` | `cb724c7122d3` | `7d32dfc88573` |
| `prl-critical-noise-research` | `3b8eaef27722` | `427b183815dc` |
| `revise-e7-manuscript` | `7b5af10a72c8` | `11d1bed2d646` |

作者、提交者、日期、提交说明保留；逐分支比较的科研文件模式、blob 与路径完全一致。因剔除了产品文件和无关祖先，提交 SHA 改变。完整映射见 [openquantum-history.json](provenance/openquantum-history.json)。历史分支应与当前工作分支一起保留；只传输当前分支不会带走其他分支独有的提交。

当前文件快照在任务分支 `c/import-research-papers` 独立提交，并建立同一提交上的 `main` 作为新仓库基线。原仓库在 `c/extract-research-papers` 独立提交迁出记录；原有科研分支不改动。

## 验证与边界

- [逐文件清单](provenance/migration-file-manifest.json)：迁移前 1,578 个文件的大小、模式、SHA256 及是否纳入当前 Git。清单是不可变的来源快照，后续科研修订不应回写它。
- [验证结果](provenance/migration-verification.json)：迁移后的内容和元数据核对、Git 历史检查、构建与只读复现摘要。
- 297 个 Python 文件通过语法解析；256 个研究 JSON 文件通过解析。未发现所检查的高置信度凭据模式；这不是完整安全审计。
- Rydberg PRA 的 89 个交付清单条目、费米 EP 的 34 个交付清单条目全部哈希匹配。
- 九份当前论文/补充材料在隔离目录编译成功，无未解析引用。八份提取文本逐字节匹配已有 PDF；SYK 仅由 `\today` 将日期从 October 27 改为 October 29，除该行外完全一致。迁移没有覆盖已有 PDF。
- 编译保留既有提示，包括 REVTeX 浮动体暂缓、nameref、部分 Underfull，以及 SYK 一处 0.91583 pt Overfull 和 PDF 书签提示。此次没有修改排版或科研正文，不把成功构建描述为零警告。
- 在新路径运行 `recovery.py --check-only`：16 组小系统 SDP 检查通过，Fourier 枚举最大误差 `1.11e-16`，CDF 最大误差 `3.33e-16`。CLARABEL 的精度提示触发既有 SCS 回退，三组涉及回退，最终状态全部为 `optimal`；没有改写研究样本。最初日志汇总因警告文本位于 JSON 前而失败，随后从已完成输出提取结果，没有重算或隐藏警告。
- 未重跑大型谱扫描、瞬子、全部科研计算或产品工程测试；迁移不构成独立科学审稿，也不改变任何项目的发表状态。

迁移快照、完整检查日志和隔离构建临时文件保存在本地 `.git/migration/`，不进入版本管理。可复核的来源清单和验证摘要已纳入 `provenance/`。未创建 GitHub 仓库、推送、投稿或公开发布。
