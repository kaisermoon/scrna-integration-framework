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

scCODA（D09_abundance 丰度分析）依赖 TF，与主环境隔离到独立 conda 环境 `scrna-sccoda`。

| 环境 | 用途 | 关键依赖 | 环境 spec |
|------|------|---------|-----------|
| `scrna-integration` | 主流水线（01-10，含 scVI/scANVI/scCRAFT） | PyTorch, scvi-tools, scanpy, rpy2 | `environment.yml` |
| `scrna-sccoda` | D09_abundance scCODA 丰度分析 | TensorFlow, scCODA | `environment-sccoda.yml`（两平台，2026-06-15 跨平台改造） |

> **跨平台改造（2026-06-15）**：`environment-sccoda.yml` 已剔除 linux 专属底层包，精简为 4 个 conda 功能包（python/pip/r-base/rpy2）+ 76 个 pip 包，与主 `environment.yml` 同风格。两平台均可直接 `conda env create -f environment-sccoda.yml`。
> **channels 说明（2026-06-17）**：channels 含 `conda-forge` + `bioconda`。两 channel 在 osx-arm64 上 R 包覆盖度不如 linux-64，但 sccoda 环境仅需 `r-base` + `rpy2`（均来自 conda-forge），bioconda 为预留（未来若需 Bioconductor R 包）。当前在两平台可正常求解。

### torch CUDA variant（预期跨平台差异，2026-06-15 登记）

| 项目 | linux-64 (CUDA) | osx-arm64 | linux-64 (CPU only) |
|------|----------------|-----------|---------------------|
| torch 来源 | pip `--index-url https://download.pytorch.org/whl/cu126` | pip PyPI | pip PyPI |
| torch 版本字符串 | `2.12.0+cu126` | `2.12.0` | `2.12.0+cpu` 或 `2.12.0` |
| env_parity 行为 | 标记为差异（预期，人工确认） | 基准（environment.yml pip PyPI） | 同 osx-arm64 基准 |

> **为什么走 pip 而不走 conda-forge**：conda-forge pytorch 2.12.0 无 cuda126 build（仅 cuda129/cuda130），
> cuda129 要求 NVIDIA driver >= R570，而项目常用 Ubuntu + A800 集群的 driver 560 仅支持到 CUDA 12.6。
> pip 提供 `torch==2.12.0+cu126` wheel，与 driver 560 兼容。
>
> **安装方式**：`conda env create -f environment.yml`（装 CPU torch），然后在 GPU 机器上运行
> `bash scripts/setup_cuda.sh` 替换为 CUDA variant。

---

### cytotrace2 独立环境（numpy 版本冲突，2026-06-16 登记）

cytotrace2-py（10d 拟时序）要求 `numpy<2.0.0`，而 scCRAFT（04 嵌入）依赖 jax 要求 `numpy>=2.0`。
两者无法共存于同一 conda 环境，仿照 scCODA 模式隔离到独立环境。

| 环境 | 用途 | 关键依赖 | 环境 spec |
|------|------|---------|-----------|
| `scrna-cytotrace2` | 10d CytoTRACE2 分化潜能 | cytotrace2-py, numpy<2.0, scanpy | `environment-cytotrace2.yml` |

> **安装方式**：`conda env create -f environment-cytotrace2.yml`
> 使用时需切换 kernel 或在 Jupyter 中选择对应 conda 环境。

### pyscenic numpy 兼容性（monkey-patch，2026-06-16 登记）

pyscenic 0.12.1 的 `transform.py:42-44` 使用 `np.object` 别名（numpy>=1.24 已移除）。
修复方式：在 `D11_grn.ipynb` 的 pyscenic 导入前执行 `np.object = object` monkey-patch，
不改动 pyscenic 源码（保护两平台兼容性）。不影响原 Mac 环境。

---

**维护规则**：本表随 env_parity 快照对比更新与 Linux 实测更新；reviewer 审 env/平台相关 PR 时核对本表与实际快照/代码一致。
