---
status: accepted
---

# 跨平台一致性：conda-lock 双平台锁文件 + 单点 OS 检测

本项目同时运行在两台机器上：**开发机 Mac（Apple Silicon，`osx-arm64`）** 与 **服务器 Linux（x86-64，`linux-64`，Alibaba Cloud Linux 3，无 GPU）**。PI 的要求是：两台机器的 conda 环境、Python/R 包版本、以及全部代码函数行为**尽可能完全一致**；确实无法兼容的，允许**最小限度**的系统检测与切换，且这种切换必须显式登记、可审计，不允许散落在各处的静默分支。

本 ADR 记录为达成此目标的环境与代码层决策，是 ADR-0007（R 桥接统一 subprocess）在多平台维度上的延伸。

## 触发的事实

1. **现有 `environment.yml` / `environment-r.yml` / `pyproject.toml` 全部用 `>=` 松约束**（`numpy>=1.24`、`scanpy>=1.10`、`r-base>=4.3` …）。松约束在「一份 spec、两台机器各自 solve」的场景下，几乎必然解出不同的版本组合——这不是装的时候小心就能避免，是 spec 本身没有锁定能力。`>=` 表达的是「库作者声明的兼容下界」（抽象依赖），不是「可复现的具体环境」（锁定依赖），两者职责不同。
2. **服务器 Linux 环境从零重建**：`miniforge3/envs/` 为空，`scrna-integration` 与 `scrna-integration-r` 都需在 Linux 上新建。这正是引入锁文件的最佳时机——重建即对齐。
3. **GPU 分歧不存在**：Linux 服务器无 NVIDIA GPU，Mac 也无 CUDA，两边 PyTorch / scVI 都跑 CPU build。跨平台最大的常见坑（CUDA vs CPU vs MPS 的不同 build）在本项目天然不存在，对齐难度大幅下降。
4. **conda 环境目录不在 Syncthing 同步范围**（项目目录在同步范围，但 `miniforge3/` 在范围外）。这是必须固守的安全状态：conda 环境内是平台相关的编译二进制，arm64 的 `.so`/`.dylib` 在 linux-64 无法运行，跨平台同步环境目录会直接损坏环境。
5. **代码层已存在不一致**：`stage2_qcd.ipynb` 从 `CONDA_PREFIX` 派生 `RSCRIPT_BIN`（带 Mac 硬编码 fallback），而 `stage7/*.ipynb` 直接写死 `RSCRIPT_BIN = "Rscript"`（靠 PATH）。同一约定两种写法，且其一含 Mac 绝对路径——换到 Linux 即失效。

## 决定

### 一、版本一致：conda-lock 双平台锁文件

1. **引入 `conda-lock` 作为唯一的版本锁定机制**。从一份多平台 spec（`conda-lock.yml` 源或现有 `environment*.yml`）为 `linux-64` 与 `osx-arm64` 两个平台**同时**求解，产出统一锁文件 `conda-{py|r}.lock.yml`。conda-lock 保证：同一个包在两个平台锁定**同一版本号**（只要该平台有对应 build），pip 依赖一并锁入（lock 覆盖 conda + pip 全栈）。
2. **职责分层**：
   - `environment.yml` / `environment-r.yml`：**人类可读的源 spec**，pin 到 Mac 当前已验证可跑的精确版本（`==`，不带 build string 以保持跨平台可移植）。这是「我们想要什么版本」的真相源。
   - `conda-{py|r}.lock.yml`：**conda-lock 生成的多平台锁文件**，含完整 transitive 依赖与每平台 build。提交进 git，是「两台机器实际装什么」的真相源。
   - `pyproject.toml`：保持 `>=` 抽象依赖（库语义，将来发 PyPI 用），**不参与环境锁定**。
3. **两台机器都从 lock 文件安装**：`conda-lock install --name scrna-integration conda-py.lock.yml`（或 `conda-lock render` 出平台专属 explicit lock 再 `conda create`）。不再用 `conda env create -f environment.yml` 直接 solve（那会重新漂移）。
4. **生成 lock 在哪台机器都行**：conda-lock 跨平台求解，不依赖运行它的机器架构。约定 **Mac 为锁文件的生成方**（已验证版本的真相源在 Mac），Linux 只消费。

### 二、对齐异常登记：无法跨平台统一的包

osx-arm64 的 bioconda 覆盖弱于 linux-64，少数 R/生信包可能在 Mac 上缺 native build。这是 PI 允许的「最小切换」唯一合法触发场景。处理规则：

1. conda-lock 求解失败的包，**优先尝试用 conda-forge/bioconda 的等价包或 noarch build** 对齐。
2. 实在无法对齐的，**在 `docs/cross-platform-exceptions.md` 显式登记**：包名、缺失平台、根因、采用的回退（如 Mac 走 Rosetta 的 `osx-64` 子环境 / 在 R 内 `BiocManager::install` / 该包仅某平台启用）、对功能的影响。
3. **登记表是异常清单，不是常态**。每新增一项 reviewer 必须质询「是否真的无法对齐」。目标是这张表尽可能短乃至为空。
4. 已知候选异常（待 Linux 重建时验证）：`monocle3` / `hdWGCNA`（本就无 conda 包，走 R 内 `BiocManager`/`remotes`，两平台同一安装路径，反而天然一致）；`r-soupx` / `bioconductor-*` 的 osx-arm64 build 可用性待实测。

### 三、代码一致：OS 检测单点收口

1. **新建 `src/scrna_integration/platform.py`**，作为全项目**唯一**的 OS/路径差异收口点。所有 `os.uname()` / `sys.platform` / `platform.system()` 判断只允许出现在此模块。
2. 首批收口对象：
   - **`rscript_bin()`**：从 `CONDA_PREFIX` 派生 `scrna-integration-r` 环境的 Rscript 路径，跨平台同一逻辑（`{envs_dir}/scrna-integration-r/bin/Rscript`，envs_dir 由 CONDA_PREFIX 上溯）。**删除所有 notebook 里的 Mac 硬编码 fallback 与写死 `"Rscript"`**，统一改调 `platform.rscript_bin()`。
   - 将来若出现真实的平台分支（临时目录约定、某包的平台条件 import 等），一律加进本模块，notebook/src 其他位置只 import 调用。
3. **判据**：notebook 与 `src` 其他文件中**不得出现** OS 判断或平台绝对路径。reviewer 审查时 grep `sys.platform` / `os.uname` / `/Users/` / `/home/` 硬编码路径，命中本模块以外的位置即 flag。

### 四、环境同步的硬约束（写入 SPEC + 项目 CLAUDE.md）

1. **conda 环境目录永不进 Syncthing / git**。只有 lock 文件、源 spec、`platform.py`、异常登记表随同步走。两台机器各自 `conda-lock install` 重建。
2. 环境隔离硬规定（ADR 既有，重申）：只装进 `scrna-integration` / `scrna-integration-r`，绝不动 base。

## Considered Options

- **手工把 `>=` 改成 `==` pin 进 yml**：部分有效，但只锁了直接依赖，transitive 依赖仍漂移；且维护靠人工 export，易遗漏。conda-lock 锁全栈且可自动重生成，是更彻底的解。保留 `==` 源 spec 作为 conda-lock 的输入，两者配合。
- **Mac 走 Rosetta 全程用 `osx-64`**：`osx-64` 与 `linux-64` 的 bioconda 覆盖最接近，对齐最容易。但 Rosetta 模拟有性能损耗（scVI 训练等重计算明显），且 PI 的 Mac 原生跑 arm64。决策为**默认 native `osx-arm64`**，仅在异常表中确有 arm64 缺包的个别工具上局部考虑 osx-64 子环境，不全程 Rosetta。
- **Docker / 单一 Linux 容器跨机统一**：能彻底消灭平台差异，但 Mac 上跑 Linux 容器同样有性能与 GPU 直通问题，且偏离「notebook 直接在本机 jupyter 跑」的项目交付形态（ADR-0009 教学透明）。容器化留作将来正式 release 时再议，当前不引入。
- **完全不锁，靠 `>=` + 人工小心**：拒绝。这正是 PI 要解决的问题——松约束无法保证一致。

## Consequences

- 新增依赖 `conda-lock`（装进 base 之外的工具环境或 pipx，不污染项目环境）；新增锁文件 `conda-py.lock.yml` / `conda-r.lock.yml` 与异常登记表 `docs/cross-platform-exceptions.md` 进 git。
- `environment.yml` / `environment-r.yml` 从 `>=` 改为 `==` 精确 pin（以 Mac 已验证版本为准）。
- README 安装说明从 `mamba env create -f environment.yml` 改为 `conda-lock install ...`；保留源 spin 说明供新增包时重生成 lock。
- 新建 `src/scrna_integration/platform.py`；现有 notebook（stage2 + stage7 全部 R-using）的 RSCRIPT_BIN 写法统一改为 `platform.rscript_bin()`。这是一次跨多 notebook 的一致性改写，走独立 PR。
- code-reviewer 新增红线：① env 相关 PR 中 `>=` 出现在源 spec（应为 `==`）或 lock 文件未同步重生成 → flag；② notebook/src 出现 OS 判断或平台绝对路径（`platform.py` 以外）→ flag；③ 新增跨平台异常未登记 `docs/cross-platform-exceptions.md` → flag。
- 全流程需在 **Linux 与 Mac 两台机器各跑一遍端到端验证**，记录版本一致性与异常表，作为本 ADR 落地的验收。
