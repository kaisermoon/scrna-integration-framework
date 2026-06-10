---
status: accepted
revised: 2026-06-10  # conda-lock 方案被 env_parity 诊断脚本方案取代（同一 ADR 内演进）
---

# 跨平台一致性：精确 pin 源 spec + env_parity 诊断脚本 + 人工对齐

本项目同时运行在两台机器上：**开发机 Mac（Apple Silicon，`osx-arm64`）** 与 **服务器 Linux（x86-64，`linux-64`，Alibaba Cloud Linux 3，无 GPU）**。PI 的要求是：两台机器的 conda 环境、Python/R 包版本、以及全部代码函数行为**尽可能完全一致**；确实无法兼容的，允许**最小限度**的系统检测与切换，且这种切换必须显式登记、可审计，不允许散落在各处的静默分支。

本 ADR 记录为达成此目标的环境与代码层决策，是 ADR-0007（R 桥接统一 subprocess）在多平台维度上的延伸。

> **修订历史**：2026-06-10，conda-lock 强锁定方案被精确 pin 源 spec + env_parity 诊断脚本 + 人工对齐方案取代（同一 ADR 内演进，不新开 ADR——原因见下）。原 conda-lock 方案记录于下文「为什么放弃 conda-lock」段与 Considered Options 中，供历史追溯。

## 触发的事实

1. **现有 `environment.yml` / `environment-r.yml` / `pyproject.toml` 全部用 `>=` 松约束**（`numpy>=1.24`、`scanpy>=1.10`、`r-base>=4.3` …）。松约束在「一份 spec、两台机器各自 solve」的场景下，几乎必然解出不同的版本组合——这不是装的时候小心就能避免，是 spec 本身没有锁定能力。`>=` 表达的是「库作者声明的兼容下界」（抽象依赖），不是「可复现的具体环境」（锁定依赖），两者职责不同。
2. **服务器 Linux 环境从零重建**：`miniforge3/envs/` 为空，`scrna-integration` 与 `scrna-integration-r` 都需在 Linux 上新建。
3. **GPU 分歧不存在**：Linux 服务器无 NVIDIA GPU，Mac 也无 CUDA，两边 PyTorch / scVI 都跑 CPU build。跨平台最大的常见坑（CUDA vs CPU vs MPS 的不同 build）在本项目天然不存在，对齐难度大幅下降。
4. **conda 环境目录不在 Syncthing 同步范围**（项目目录在同步范围，但 `miniforge3/` 在范围外）。这是必须固守的安全状态：conda 环境内是平台相关的编译二进制，arm64 的 `.so`/`.dylib` 在 linux-64 无法运行，跨平台同步环境目录会直接损坏环境。
5. **代码层已存在不一致**：`stage2_qcd.ipynb` 从 `CONDA_PREFIX` 派生 `RSCRIPT_BIN`（带 Mac 硬编码 fallback），而 `stage7/*.ipynb` 直接写死 `RSCRIPT_BIN = "Rscript"`（靠 PATH）。同一约定两种写法，且其一含 Mac 绝对路径——换到 Linux 即失效。

## 决定

### 一、版本一致：精确 pin 源 spec + env_parity 诊断脚本 + 人工对齐

1. **源 spec 用精确 `==` pin**（`environment.yml` / `environment-r.yml`）。pin 到 Mac 当前已验证可跑的精确版本（不带 build string 以保持跨平台可移植）。这是「期望基准」——两机装出来应尽可能一致。`pyproject.toml` 保持 `>=` 抽象依赖（库语义，不参与环境锁定）。
2. **两机各自 `conda env create -f environment.yml`（或 mamba）直接安装**。精确 pin 已尽量保证一致。
3. **诊断工具 `scripts/env_parity.py`**（纯 Python 标准库 + subprocess 调 conda）：
   - `snapshot`：感知当前机器（通过 `platform_tag()` 识别 linux-64/osx-arm64）、导出两个 conda 环境的完整包清单（`conda list --json`）、写入 `docs/env-snapshots/{platform_tag}.json`。
   - `compare`：对比两份快照（默认 linux-64 vs osx-arm64），输出版本不一致、仅 A 有、仅 B 有的包三类差异表。完全一致则打印 ✓ 并返回 0；有差异返回 1（便于 CI/hook 感知）。
   - **硬约束：脚本只诊断报告差异，绝不自动改环境**（对齐决策归人，不外包给脚本）。
4. **人工对齐流程**：`compare` 出差异后，由人在落后那台执行 `conda install pkg=version` 或调 environment.yml 后重建。对齐后重新 `snapshot` 更新快照，提交 git 供下次对比。

### 二、对齐异常登记：无法跨平台统一的包

osx-arm64 的 bioconda 覆盖弱于 linux-64，少数 R/生信包可能在 Mac 上缺 native build。这是 PI 允许的「最小切换」唯一合法触发场景。处理规则：

1. 优先尝试用 conda-forge/bioconda 的等价包或 noarch build 对齐。
2. 实在无法对齐的，**在 `docs/cross-platform-exceptions.md` 显式登记**：包名、缺失平台、根因、采用的回退（如 Mac 走 Rosetta 的 `osx-64` 子环境 / 在 R 内 `BiocManager::install` / 该包仅某平台启用）、对功能的影响。
3. **登记表是异常清单，不是常态**。每新增一项 reviewer 必须质询「是否真的无法对齐」。目标是这张表尽可能短乃至为空。
4. 已知候选异常（待 Linux 重建时验证）：`monocle3` / `hdWGCNA`（本就无 conda 包，走 R 内 `BiocManager`/`remotes`，两平台同一安装路径，反而天然一致）；`r-soupx` / `bioconductor-*` 的 osx-arm64 build 可用性待实测。

### 三、代码一致：OS 检测单点收口

1. **`src/scrna_integration/platform.py`**，作为全项目**唯一**的 OS/路径差异收口点。所有 `os.uname()` / `sys.platform` / `platform.system()` 判断只允许出现在此模块。
2. 首批收口对象：
   - **`rscript_bin()`**：从 `CONDA_PREFIX` 派生 `scrna-integration-r` 环境的 Rscript 路径，跨平台同一逻辑（`{envs_dir}/scrna-integration-r/bin/Rscript`，envs_dir 由 CONDA_PREFIX 上溯）。**删除所有 notebook 里的 Mac 硬编码 fallback 与写死 `"Rscript"`**，统一改调 `platform.rscript_bin()`。
   - **`platform_tag()`**：返回标准化平台标识（`linux-64` / `osx-arm64` / `osx-64`），供快照文件命名与跨机器对比。
   - 将来若出现真实的平台分支（临时目录约定、某包的平台条件 import 等），一律加进本模块，notebook/src 其他位置只 import 调用。
3. **判据**：notebook 与 `src` 其他文件中**不得出现** OS 判断或平台绝对路径。reviewer 审查时 grep `sys.platform` / `os.uname` / `/Users/` / `/home/` 硬编码路径，命中本模块以外的位置即 flag。

### 四、环境同步的硬约束（写入 SPEC + 项目 CLAUDE.md）

1. **conda 环境目录永不进 Syncthing / git**。只有源 spec、`platform.py`、异常登记表、快照 JSON 随同步走。两台机器各自 `conda env create` 安装。
2. 环境隔离硬规定（ADR 既有，重申）：只装进 `scrna-integration` / `scrna-integration-r`，绝不动 base。

## 为什么放弃 conda-lock

最初方案（2026-06-08）采用 `conda-lock` 为双平台生成统一锁文件，两机从 lock 安装以保证版本完全一致。PI 在 2026-06-10 决定放弃此方向，理由如下：

- **维护成本高**：项目会持续装新包（下游模块不断加入新 R/Bioc 工具、新 Python 包），每次都要重生成 lock 文件并两机重装。锁文件机制适合「依赖冻结后长期不变」的场景（如论文 camera-ready 版本归档），不适合单 PI 持续演进的科研环境。
- **僵硬**：lock 文件锁死所有传递依赖的 build string，即使在 conda-forge 上新 build 修复了 bug，也需要人工介入才能更新。精确 pin 源 spec 保留了一定的 conda 求解灵活性——版本号保一致，但 conda 可以在同版本的不同 build 之间选择最适配当前平台的。
- **诊断 + 人工对齐更适合本项目规模**：两台机器、两个环境、一位 PI。差异量小（精确 pin 已消除了绝大部分差异），`env_parity.py compare` 出的差异人眼扫一遍就能决定。不需要为一两个包的微小差异引入一整套 lock 重生成流水线。

## Considered Options

- **conda-lock 双平台锁文件**（原方案，2026-06-08）**→ 2026-06-10 被取代**：能彻底锁死全栈版本，但维护成本高、不适合持续演进。保留此记录供将来（如论文提交前需要完全冻结时）重新评估。
- **手工把 `>=` 改成 `==` pin 进 yml**：采纳。这是当前方案的基石——源 spec 用 `==` 精确 pin，传递依赖虽有微小漂移空间但已在 conda 求解器的约束范围内。配合 `env_parity compare` 诊断残余差异，人决定是否对齐。
- **Mac 走 Rosetta 全程用 `osx-64`**：`osx-64` 与 `linux-64` 的 bioconda 覆盖最接近，对齐最容易。但 Rosetta 模拟有性能损耗（scVI 训练等重计算明显），且 PI 的 Mac 原生跑 arm64。决策为**默认 native `osx-arm64`**，仅在异常表中确有 arm64 缺包的个别工具上局部考虑 osx-64 子环境，不全程 Rosetta。
- **Docker / 单一 Linux 容器跨机统一**：能彻底消灭平台差异，但 Mac 上跑 Linux 容器同样有性能与 GPU 直通问题，且偏离「notebook 直接在本机 jupyter 跑」的项目交付形态（ADR-0009 教学透明）。容器化留作将来正式 release 时再议，当前不引入。
- **完全不锁，靠 `>=` + 人工小心**：拒绝。这正是 PI 要解决的问题——松约束无法保证一致。

## Consequences

- 源 spec（`environment.yml` / `environment-r.yml`）从 `>=` 改为 `==` 精确 pin（以 Mac 已验证版本为准），作为期望基准。
- 新建 `scripts/env_parity.py`（snapshot + compare 两子命令）；快照目录 `docs/env-snapshots/`（每台机器一份 JSON，进 git）。
- 新建 `src/scrna_integration/platform.py`（含 `rscript_bin()` + `platform_tag()`）；现有 notebook（stage2 + stage7 全部 R-using）的 RSCRIPT_BIN 写法统一改为 `platform.rscript_bin()`。这是一次跨多 notebook 的一致性改写，走独立 PR。
- 异常登记表 `docs/cross-platform-exceptions.md` 保持，维护规则从"随 conda-lock 重生成"改为"随 env_parity 快照对比更新"。
- README 安装说明从 `conda-lock install ...` 改为 `conda env create -f environment.yml`。
- code-reviewer 红线更新：① env 相关 PR 中 `>=` 出现在源 spec（应为 `==`）→ flag；② notebook/src 出现 OS 判断或平台绝对路径（`platform.py` 以外）→ flag；③ 新增跨平台异常未登记 `docs/cross-platform-exceptions.md` → flag；④ 快照 JSON 中包名/版本与源 spec 不一致且 PR 未说明 → flag。
- 全流程需在 **Linux 与 Mac 两台机器各跑一遍端到端验证**，记录版本一致性与异常表，作为本 ADR 落地的验收。
