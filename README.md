# CASCAT

**CASCAT — Coverage-Aware Search for Catalytic Adsorbate Topologies**，用于催化金属表面在不同吸附数量下的构型搜索。用户提供 clean 金属表面的 POSCAR / CONTCAR 和吸附分子，CASCAT 在本机完成逐覆盖度搜索。项目基于已运行的 ACAT + MACE 流程，提供独立 Python 包与 `cascat` 命令；不需要原项目目录、个人 conda 路径或 VASP 安装。

## 安装

### 让助手完成安装与运行

本仓库附带 [run-cascat-search skill](skill/run-cascat-search/SKILL.md)。下载或克隆完整仓库后，在助手中打开本项目目录，使用：

```text
请读取 skill/run-cascat-search/SKILL.md，按其中流程帮我安装 CASCAT 并运行：
表面文件是 /path/to/POSCAR，吸附 CO，搜索到 N=4。
每个覆盖度从前一覆盖度的 2 个构型出发；
单点后取前 80%，粗优化后取前 50%，精修后保存前 10 个结果。
```

skill 会指导助手检查环境、安装依赖、核对输入、执行搜索、跟进运行并检查最终文件。缺少必要的分子/数量/输入文件时会先询问；已有信息不会重复询问。也可仅要求安装、只生成候选、续跑或检查结果。

skill 放在可见的 `skill/` 目录中；使用时明确让助手读取上述 `SKILL.md`，不依赖自动发现。只执行 `pip install .` 不会将 skill 安装到用户全局目录。

上传 GitHub 时请一并上传 `skill/` 目录和 `MANIFEST.in`；skill 的说明、安装参考和检查脚本都在其中。源代码发行包也会包含这些文件。

### 手动安装

使用 Python 3.11 或 3.12。在本项目目录执行（macOS / Linux；Windows 推荐 WSL2）：

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
cascat molecules
```

首次安装依赖和首次运行 MACE 需要联网下载。模型使用 MACE-MH-1 的 `oc20_usemppbe` head；缓存后可离线运行。默认 CPU、1 个优化进程、每进程 2 线程。每个进程单独加载模型；内存足够时可用 `--workers 4`。

## 最简单的使用方法

```bash
cascat run /path/to/POSCAR --ads CO --n 4 --parents 2 \
  --screen-keep 80% --relax-keep 50% --final-keep 10 --output results/CO
```

这会从 clean slab 开始依次搜索 1、2、3、4 个 CO。`--n` 是目标分子数；`--ads` 是吸附物名称。**父构型数量与各阶段筛选数量分别控制**：

| 参数 | 含义 | 接受的值 | 默认 |
| --- | --- | --- | --- |
| `--parents` | N=n 从 N=n−1 的最终保留结果中，取最低能的多少个构型出发 | 正整数或 `all` | `1` |
| `--screen-keep` | 初始候选去重、单点能量排序后，多少个进入粗优化 | 正整数、百分比或 `all` | `80%` |
| `--relax-keep` | 粗优化排序去重后，多少个进入精修 | 正整数、百分比或 `all` | `50%` |
| `--final-keep` | 精修排序去重后，保存多少个合格结果，其中最低能的是 BEST | 正整数、百分比或 `all` | `all` |

上例中，每个 N 都先保留单点排名前 80% 进入粗优化，再从粗优化后的唯一有效构型中保留前 50% 进入精修，最终保存最低能的至多 10 个合格结果。下一覆盖度只使用这些结果中的最低能至多 2 个父构型，各增加一个分子继续搜索。N=1 从 clean slab 出发。

`--parents` 只设一个值，所有覆盖度共用，不影响三个阶段的筛选比例或数量。若最终只保存了 3 个结果，即便 `--parents 5` 也只能使用这 3 个。希望从更多结果出发，需要同时让 `--final-keep` 保留足够的结果。

各阶段也可使用不同的固定数量：

```bash
cascat run CONTCAR --ads OH --n 4 --parents 2 \
  --screen-keep 1024 --relax-keep 256 --final-keep 10 --output results/OH
```

百分比以**当前阶段去重后的合格结构数**为分母，各阶段独立计算，向上取整。例如 11 个取 80% 得到 9 个，7 个取 50% 得到 4 个。有候选时至少保留 1 个；没有合格结构时不会凭空补足。整数从 1 起，不足指定数量时使用实际数量。百分比范围为大于 0、至多 100%，支持小数。

各参数都可单独设为 `all`，表示该处不限制数量。三个阶段仍执行去重、完整性判断和最终收敛判断，不额外套用能量窗口或隐藏数量上限。`--parents all` 表示使用上一覆盖度**最终保留的全部结果**；要让所有精修合格结果都成为父构型，同时设 `--final-keep all`。较大覆盖度的计算量可能迅速增长。

0.3.0 起移除含义混合的 `--keep`，改用上述四个独立参数。旧版结果请保留，使用新的输出目录运行新规则。

默认原样使用输入晶胞，不自动扩成 2×2。需要扩胞时指定 `--repeat 2 2 1`。表面须位于 xy 平面，z 方向有真空，吸附在上表面。输入应为不含已有吸附物的 clean 金属 slab；基底所有原子固定，仅吸附物参与优化。

默认 `--surface auto` 使用原工作流的 AutoAdsorbate / ACAT CustomSurface 路径识别表面。已知 FCC 晶面可用 `--surface 111` 等标签；不要给 BCC/HCP 输入指定 FCC 标签。支持范围仍受位点识别算法限制，不承诺任意输入都能正确识别。

## 可直接运行的例子

`examples/POSCAR` 为 3×3×4 层 Cu(111) 示例，也可重新生成：

```bash
python examples/make_slab.py
cascat run examples/POSCAR --ads H --n 2 --parents 1 \
  --screen-keep 80% --relax-keep 50% --final-keep all --surface 111 --output results/demo
```

快速试跑可将 `--screen-keep` 和 `--relax-keep` 都设为 `1`。这会缩小搜索空间。

只检查候选生成、暂不加载 MACE：

```bash
cascat run examples/POSCAR --ads CO2 --n 4 --parents 2 --surface 111 --output results/CO2 --generate-only
```

该模式只生成**下一未完成覆盖度**的全部去重候选，暂不按 `--screen-keep` 截断，因为排序需要 MACE 单点能量。后续覆盖度需要上一轮优化后的父构型。去掉 `--generate-only` 后执行同一命令即可正式运行。

## 分子库

用 `cascat molecules` 列出模板与内部命名。当前包含 28 种，覆盖 HER、OER 及常见 CO₂RR C1 路径和部分 C2 物种：

| 范围 | 模板 |
| --- | --- |
| HER 与水相关 | H、H2、H2O、OH |
| OER | O、OH、OOH、O2、H2O、H2O2 |
| CO₂RR C1 | CO2、COOH、HCOO、HCOOH、CO、CHO、COH、CHOH、CH2O、CH2OH、CH3O、CH3OH、C、CH、CH2、CH3、CH4 |
| CO₂RR 部分 C2 | CH2CO、CH3CO、CH3COOH |

这里的模板是放置在表面的中性几何片段，不赋予溶液物种的电荷态。HCOO 是 formate 几何模板；H⁺、OH⁻、电子和电解液不作为独立的吸附搜索对象。一次任务使用一种吸附物，不做混合吸附。

CO₂RR 的完整物种集合取决于研究路径。本版不声称包含“全部中间体”：例如完整 C–C 偶联网络、乙烯/乙醇等路径还需扩展模板并验证。新增条目在 `src/cascat/molecules.py` 中集中维护；新增模板必须通过生成与识别测试，不能只添加化学式名称。

模板复用 ACAT 几何与结合顺序命名，例如 H2O → OH2、CH3O → OCH3。多齿分子按原子集合统计分子数，避免把多个占据位点算作多个分子。识别使用 ACAT 的连接与位点判断；模板能生成不表示它在特定表面优化后一定稳定。参见 [ACAT 构建文档](https://asm-dtu.gitlab.io/acat/build.html) 和 [覆盖度识别文档](https://asm-dtu.gitlab.io/acat/modules.html)。

## 搜索与输出

沿用原脚本的主要计算步骤：ACAT 增量生成与对称去重 → float32 单点能量与排序选择 → 粗优化 → 吸附图去重与排序选择 → float64 精修 → 最终去重与结果选择 → 独立选择下一覆盖度的父构型。短优化与粗优化指同一阶段，完整优化与精修指同一阶段。

粗优化默认 10 步，精修累计最多 120 步，收敛阈值 0.05 eV/Å。各阶段分别使用对应的筛选参数。H/O 保留原有 H₂/O₂ 形成排除阈值。运行参数用 `cascat run --help` 查看。

每轮输出位于 `results/任务名/N_1/`、`N_2/` 等目录：

- `BEST.vasp`：该轮最低合格构型。
- `final_results.traj` / `.xyz`：按精修能量排序、由 `--final-keep` 选中的结果，第一帧对应 BEST。
- `growth_parents.traj` / `.xyz`：由 `--parents` 从最终结果中选出的下一覆盖度父构型。
- `candidate_pool.traj` / `.xyz`：初始候选。
- `ranking.csv`：完整排名、收敛与图去重信息；`retained_final` 和 `growth_parent` 分别标记最终保留结果与父构型。未选中的结构仍留在阶段轨迹和排名中用于审计。
- `pre_relaxation_ranking.csv` 与各阶段轨迹：中间结果。
- `metadata.json`：状态、实际保留数量和筛选记录。

任务根目录保存 `slab.vasp` 和 `search.json`（输入文件哈希、参数与依赖版本）。重跑相同命令会复用状态为 `complete` 的轮次；未完成轮次重新计算。修改分子、输入、保留数、总轮数或计算参数，请使用新的输出目录。程序不自动覆盖不同配置的结果。

同一输出目录禁止多个任务同时运行。正常退出会移除 `.running.lock`；强制杀死进程后，确认旧进程已停止再手动移除该锁文件。若候选无法继续增长，程序明确报错；本版不自动改用未选中的旧候选，以保证用户指定的父构型范围。

## 验证与限制

```bash
python -m unittest discover -s tests -v
```

测试覆盖全部库条目的单分子候选生成与完整识别、超过两个父构型的保留规则，以及输入变化时拒绝误续跑。实际 MACE 两轮验证与依赖信息见 `VALIDATION.md`。

本版优先面向原工作流的金属表面。含 O/C/H/N/S/F 等基底元素的氧化物、碳材料等不支持，程序会拒绝以免与吸附物混淆。未在所有操作系统、新环境或所有分子上完成生产级搜索验证。

这是有限候选与有限父构型的构象搜索，不保证全局最低能。MACE 能量只用于相同吸附物、相同数量体系内部排序；不代表恒电势、溶剂化自由能或最终 DFT 结果。原有 Cu/Pt 数据与脚本仍保留在上级目录；本项目可单独复制分发，无需带上历史计算数据。
