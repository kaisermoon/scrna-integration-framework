# 跨平台对齐异常登记表

> **ADR-0010**。本表登记**无法在 `linux-64` 与 `osx-arm64` 之间用同一来源/版本对齐**的包，及其回退方案。
>
> **这是异常清单，不是常态。** 每新增一项，code-reviewer 必须质询「是否真的无法对齐」。目标：本表尽可能短乃至为空。两台机器的 conda 环境、包版本、代码函数行为应尽可能完全一致；只有确实无法跨平台用同一来源/版本对齐时，才允许在此登记最小切换。

## 一、无法跨平台用同一来源/版本对齐的包

| 包 | 缺失平台 | 根因 | 回退方案 | 功能影响 | 登记日期 |
|---|---|---|---|---|---|
| `r-wgcna` | osx-arm64 | bioconda 仅发布 linux-64 build，osx-arm64 无包；放进 environment-r.yml conda 层会让 Mac 求解失败 | 不进 conda，改在 R 内 `install.packages("WGCNA")`（CRAN 两平台都有，源码编译，依赖 `compilers` 元包）。归入第二类统一管理 | stage7 hdWGCNA 依赖 WGCNA；R 内装后两平台一致 | 2026-06-10 |

> **本机 Linux 当前偏差**：2026-06-10 首次重建时 operator 用 `bioconda r-wgcna=1.74` 装进了 `scrna-integration-r`（因 linux-64 有包）。Mac 无此 build，故规范改为两平台统一走 CRAN R 内装。本机 Linux 的 conda r-wgcna 1.74 暂留可用，**PR-X3 双机对齐时改为 R 内装以与 Mac 一致**（或确认 CRAN 版本与 1.74 一致即不动）。

## 二、不进 conda、两平台同路径安装的包（非异常，仅记录）

这些包本就无 conda 包（或 conda 仅单平台），重建后在 R 内统一安装，两平台路径一致，列此为安装可复现：

| 包 | 安装方式 | 用途 | 本机 Linux 实装版本 |
|---|---|---|---|
| `WGCNA` | `install.packages("WGCNA")`（CRAN） | hdWGCNA 依赖 + stage7 共表达 | 当前为 conda 1.74（见第一类，待改 CRAN） |
| `monocle3` | `BiocManager::install("monocle3")` 或 `remotes::install_github("cole-trapnell-lab/monocle3")` | stage7 拟时序轨迹（subprocess Rscript） | v1.4.27（依赖 r-leidenbase 已 conda 预备） |
| `hdWGCNA` | `remotes::install_github("smorabit/hdWGCNA")` | stage7 共表达模块（subprocess Rscript） | v0.4.11（dev 分支，curl 下 tarball 装） |
| `CellChat` | `remotes::install_github("jinworks/CellChat")` | stage7 细胞通讯（subprocess Rscript） | v2.2.0.9001 |
| `scCRAFT` | `git clone https://github.com/ch2343/scCRAFT && cd scCRAFT && pip install .`（非 PyPI） | 04_embedded 可选嵌入方法（anchor-free VAE 批次校正） | 2026-06-13 | v0.1.0（git clone） |

> 这些 R 包两平台都从源码编译，依赖 environment-r.yml 的 `compilers` 元包提供工具链。本机 Linux 装它们额外补了系统库：hdf5 / cairo / udunits2 / gdal / proj / geos + bioconductor-rhdf5\* + CRAN ggraph / tidygraph / enrichR / tester。安装命令与版本固定后写入 `scripts/install_r_github_pkgs.R`（PR-X3）供两机一键复现。
> scCRAFT 装后需确认 TF 未被其依赖拉回（`pip uninstall -y tensorflow tensorflow-probability keras 2>/dev/null || true`）。两平台均需执行此检查（详见 ADR-0012）。

## 二点五、in-process R 桥接的双 R 版本并存（预期设计，非异常）

| 环境 | R 版本 | 用途 |
|---|---|---|
| `scrna-integration`（Python） | 4.5.x | conda 装 rpy2=3.6.7 自动拉入匹配 r-base，专供 rpy2/anndata2ri 的 in-process Python↔R 桥接 |
| `scrna-integration-r`（R） | 4.4.3 | subprocess Rscript 调用（SoupX/DESeq2/monocle3/CellChat/hdWGCNA 等重型 R） |

> 两 R 版本**职责分离、各自独立**——桥接走 Python 环境内 4.5.x，subprocess 走 R 环境 4.4.3，是设计而非冲突。rpy2 走 conda 不走 pip 的原因见 environment.yml 头部（Linux pip rpy2 链接 R 库失败）。rpy2 3.6.7 在两平台 conda-forge 都有同版本。
>
> **本机 Linux 当前偏差**：首次重建时 pip 残留了 `rpy2-rinterface`(3.6.6)/`rpy2-robjects`(3.6.5) 子包与 conda rpy2 主体混存。**重建/PR-X3 应确保纯 conda 装 rpy2**，不让 pip 再装子包。anndata2ri==2.0 仍走 pip（conda-forge 无此包）。

## 三、代码层平台分支（收口在 platform.py）

所有 OS 检测收口到 `src/scrna_integration/platform.py`（ADR-0010）。当前收口项：

| 函数 | 平台差异 | 处理 |
|---|---|---|
| `rscript_bin()` | Rscript 可执行路径 | 从 `CONDA_PREFIX` 上溯派生 `{envs_dir}/scrna-integration-r/bin/Rscript`，两平台同一逻辑，无硬编码绝对路径 |

> notebook 与 src 其他位置**禁止**出现 `sys.platform` / `os.uname` / `/Users/` / `/home/` 等平台分支或绝对路径。新增平台差异一律加进 platform.py 并在此登记。
| `platform.env_check()` | 主环境误装 TF 检测 | 从 `CONDA_PREFIX` 识别当前环境 + `importlib.util.find_spec` 检测 tensorflow/keras，只诊断不修环境（ADR-0010 延伸） | 2026-06-13 |
| `platform.detect_device()` | 计算设备（CUDA / MPS / CPU）三环境自适应 | 三环境自适应：auto 模式 CUDA→gpu > Mac scVI/scANVI→cpu > cpu；显式 `DEVICE=cuda / mps / cpu` 覆盖；scVI/scANVI Mac 默认 CPU（MPS 数值稳定性未验证），scCRAFT 恒 CPU；见 ADR-0013 | 2026-06-14 |

## 四、环境级隔离与约束（两平台统一）

### 主环境禁装 TensorFlow（ADR-0012）

TF 2.21.0 的 oneDNN/MKL runtime 与 PyTorch 2.12.0 的 MKL 在 backward pass 冲突，导致 scCRAFT/scVI 训练 C++ 层 segfault（exit 139）。两平台主环境 `scrna-integration` 均不得安装 `tensorflow`、`tensorflow-probability`、`keras`、`tensorboard`。`platform.env_check()` 在每次 notebook setup cell 自动检测此项。

### scCODA 独立环境 `scrna-sccoda`

scCODA（11_abundance 丰度分析）依赖 TF，与主环境隔离到独立 conda 环境 `scrna-sccoda`。

| 环境 | 用途 | 关键依赖 | 环境 spec |
|------|------|---------|-----------|
| `scrna-integration` | 主流水线（01-10，含 scVI/scANVI/scCRAFT） | PyTorch, scvi-tools, scanpy, rpy2 | `environment.yml` |
| `scrna-sccoda` | 11_abundance scCODA 丰度分析 | TensorFlow, scCODA | `environment-sccoda.yml`（当前 linux-64 专属） |

> **跨平台差异**：`environment-sccoda.yml` 当前为 linux-64 专属（pin 了 `__linux` virtual package 及 linux 底层库如 `libgcc-ng`、`ld_impl_linux-64`）。Mac 需手动创建：`conda create -n scrna-sccoda python=3.11 -y && pip install sccoda`。详见 `docs/MAC-SYNC.md`。

---

**维护规则**：本表随 env_parity 快照对比更新与 Linux 实测更新；reviewer 审 env/平台相关 PR 时核对本表与实际快照/代码一致。
