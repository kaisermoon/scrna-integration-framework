# 跨平台对齐异常登记表

> **ADR-0010**。本表登记**无法在 `linux-64` 与 `osx-arm64` 之间用 conda-lock 统一锁定同一版本**的包，及其回退方案。
>
> **这是异常清单，不是常态。** 每新增一项，code-reviewer 必须质询「是否真的无法对齐」。目标：本表尽可能短乃至为空。两台机器的 conda 环境、包版本、代码函数行为应尽可能完全一致；只有 conda-lock 确实解不出跨平台同版本时，才允许在此登记最小切换。

## 一、conda-lock 无法跨平台对齐的包

| 包 | 缺失平台 | 根因 | 回退方案 | 功能影响 | 登记日期 |
|---|---|---|---|---|---|
| _（待 PR-X3 在 Linux 重建时实测填入）_ | | | | | |

## 二、不进 conda、两平台同路径安装的包（非异常，仅记录）

这些包本就无 conda 包，重建后在 R 内统一安装，两平台路径一致，**不算跨平台异常**，列此仅为安装可复现：

| 包 | 安装方式 | 用途 | 备注 |
|---|---|---|---|
| `monocle3` | `BiocManager::install("monocle3")` 或 `remotes::install_github("cole-trapnell-lab/monocle3")` | stage7 拟时序轨迹（subprocess Rscript） | 依赖 r-leidenbase（已在 environment-r.yml 预备） |
| `hdWGCNA` | `remotes::install_github("smorabit/hdWGCNA")` | stage7 共表达模块（subprocess Rscript） | |
| `CellChat` | `remotes::install_github("jinworks/CellChat")` | stage7 细胞通讯（subprocess Rscript） | |

> 这三个 R 包在两平台都从源码编译安装，依赖 environment-r.yml 里的 `compilers` 元包提供工具链。安装命令与版本固定后写入 `scripts/install_r_github_pkgs.R`（PR-X3）供两机一键复现。

## 三、代码层平台分支（收口在 platform.py）

所有 OS 检测收口到 `src/scrna_integration/platform.py`（ADR-0010）。当前收口项：

| 函数 | 平台差异 | 处理 |
|---|---|---|
| `rscript_bin()` | Rscript 可执行路径 | 从 `CONDA_PREFIX` 上溯派生 `{envs_dir}/scrna-integration-r/bin/Rscript`，两平台同一逻辑，无硬编码绝对路径 |

> notebook 与 src 其他位置**禁止**出现 `sys.platform` / `os.uname` / `/Users/` / `/home/` 等平台分支或绝对路径。新增平台差异一律加进 platform.py 并在此登记。

---

**维护规则**：本表随 conda-lock 重生成与 Linux 实测更新；reviewer 审 env/平台相关 PR 时核对本表与实际 lock/代码一致。
