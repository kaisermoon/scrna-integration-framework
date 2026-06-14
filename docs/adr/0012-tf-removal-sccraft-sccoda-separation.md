# ADR-0012: 主环境移除 TensorFlow + scCRAFT 集成 + scCODA 环境隔离

- **Status**: Accepted
- **Date**: 2026-06-13
- **Supersedes**: None
- **Related**: ADR-0010 (跨平台一致性), ADR-0008 (吸收 student-code 并重写)

## Context

在集成 scCRAFT（anchor-free VAE 批次校正方法，作为 04_embedded 可选嵌入方法）时，训练阶段出现 segfault（exit 139），scCRAFT 代码本身未做任何修改。

根因排查过程：

1. 初始假设：scCRAFT 依赖 jax，可能与 PyTorch 冲突。但测试排除了此假设——单独 patch scCRAFT 的 `import jax` 不生效，segfault 依旧。
2. 最终根因（实测确认）：`opt_einsum` 被 torch/scvi 导入链触发后会 lazy-import tensorflow backend。环境中的 TensorFlow 2.21.0 初始化 oneDNN/MKL runtime，与 PyTorch 2.12.0 的 MKL 在 backward pass 产生冲突，导致 C++ 层 segfault。
3. 验证证据：移除 TensorFlow 后，**原始未改动的** scCRAFT 训练正常完成（30 epochs，loss 72.9 → 19.0，`X_scCRAFT` 无 NaN）。

同时，scCODA（11_abundance 丰度分析依赖的细胞组成差异分析工具）依赖 TensorFlow，与 scCRAFT/PyTorch 不能在同一 conda 环境共存。

## Decision

### 1. 主环境 `scrna-integration` 移除 TensorFlow 全家桶

从 `environment.yml` 移除：`tensorflow`、`tensorflow-probability`、`keras`、`tensorboard`。

### 2. scCRAFT 通过 git clone + pip install 装入主环境

- scCRAFT 不是 PyPI 包，通过 `git clone + pip install .` 安装到主环境
- 装后必须确认 TF 未被 scCRAFT 的依赖拉回（`pip uninstall -y tensorflow tensorflow-probability keras 2>/dev/null || true`）
- 装 scCRAFT 后重新生成环境快照（`env_parity.py snapshot`），供跨机对比

### 3. scCODA 隔离到独立环境 `scrna-sccoda`

- 新建 `scrna-sccoda` conda 环境，含 TensorFlow + scCODA
- 通过 `environment-sccoda.yml` 维护其依赖 spec
- 主环境 `environment.yml` 不再包含 scCODA 或 TF 相关包
- 仅在需要运行 11_abundance scCODA 分析时切换到该环境

### 4. 新增 `platform.env_check()` 自动诊断

在 `platform.py` 新增 `env_check()` 函数。每次运行主流水线 notebook（03-06）的 setup cell 自动调用，检测：

- 当前平台（linux-64 / osx-arm64 / osx-64）
- 当前 conda 环境名称（是否为 `scrna-integration`）
- 主环境中是否误装了 TensorFlow（打印调整命令）

遵循 ADR-0010 哲学：**只诊断问题、给出建议、不自动改环境**。环境调整由人工确认后执行。

## Consequences

### 正面

- scCRAFT 和 scVI/scANVI 在主环境稳定可用，不再发生训练 segfault
- `env_check()` 让环境问题在运行前暴露，避免训练到一半才 crash
- scCODA 环境隔离避免了 TF/PyTorch 冲突，且职责边界清晰

### 负面

- scCODA 需要切换到独立环境才能运行（增加了操作步骤）
- `environment.yml` 与 `environment-sccoda.yml` 两套 spec 需各自维护
- `environment-sccoda.yml` 当前为 linux-64 专属（含 linux 底层库 pin），Mac 需手动创建（见 `docs/MAC-SYNC.md`）

### 跨平台影响

- **Mac 同步要求**：主环境同样需要移除 TF（见 `docs/MAC-SYNC.md`）
- **scCRAFT on Mac**：Apple Silicon 上的 GPU 加速依赖（torch MPS / jax-metal）与 Linux CUDA 不同，需按 Mac 适配
- **scCODA on Mac**：Mac 上不能直接 `conda env create -f environment-sccoda.yml`（含 linux-64 专属包），需手动 `conda create + pip install sccoda`

### 与 ADR-0010 的关系

`env_check()` 是 ADR-0010"代码层 OS 检测收口到 `platform.py`"的延伸——新增环境自检能力，遵循同样的"只诊断不自动改环境"原则。所有环境诊断逻辑收口在单一模块，notebook 与 src 其他位置不出现重复的环境检测代码。

## Alternatives Considered

1. **保留 TF，patch scCRAFT 绕过冲突**：排除了。scCRAFT 的 segfault 发生在 PyTorch backward pass 的 C++ 层，非 Python 层面可控制；且 TF 初始化发生在 `opt_einsum` 的 lazy import 阶段（比 scCRAFT 代码更早），patch scCRAFT 无法阻止。
2. **用环境变量禁用 oneDNN/MKL**：未采纳。`TF_ENABLE_ONEDNN_OPTS=0` 等环境变量可以部分抑制 TF 的 oneDNN 初始化，但不能完全消除冲突，且可能影响其他依赖 MKL 的包（如 numpy/scipy）的性能。
3. **全程保留 TF，scCRAFT 用单独环境**：未采纳。scCRAFT 是 04_embedded 的核心替代嵌入方法，需与 scVI/scANVI 在同一个 notebook 中对比，同环境才能保证参数（epochs/batch/latent dim）完全一致。scCODA 仅在 11_abundance 使用，且与 scVI 下游无紧密的参数交叉，隔离影响最小。
