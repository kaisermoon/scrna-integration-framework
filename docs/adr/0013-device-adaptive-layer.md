---
status: accepted
---

# ADR-0013: Device 自适应层（Mac-MPS / Linux-CPU / Ubuntu-CUDA 三环境）

- **Status**: Accepted
- **Date**: 2026-06-14
- **Supersedes**: None
- **Related**: ADR-0010 (跨平台一致性), ADR-0012 (主环境移除 TensorFlow)

## Context

ADR-0010 曾判断"GPU 分歧不存在"——当时项目仅跑在两台无 GPU 的机器上（Mac MPS 未使用 + Linux 无 NVIDIA GPU），PyTorch / scVI 在所有环境均走 CPU，设备差异不是一个需要处理的问题。

现在引入第三台机器——**Ubuntu + CUDA GPU 生产环境**。同一套代码需要在三种设备条件下运行：

| 环境 | 操作系统 | 可用设备 | scVI train accelerator |
|------|----------|----------|------------------------|
| Mac 开发机 | macOS (osx-arm64) | MPS backend (Apple Silicon) | `"mps"` 或 `"cpu"` |
| Linux 服务器 | Alibaba Cloud Linux 3 (linux-64) | 无 GPU | `"cpu"` |
| Ubuntu 生产机 | Ubuntu (linux-64) | NVIDIA CUDA GPU | `"gpu"`（scvi-tools 合法值，非 `"cuda"`） |

此外，三个嵌入方法对设备的处理各有约束：

1. **scVI / scANVI (scvi-tools 1.4.2)**：`model.train()` 通过 `accelerator` + `devices` 参数控制设备（旧参数 `use_gpu` 已移除）。合法 `accelerator` 值：`"cpu"` / `"gpu"` / `"mps"` / `"auto"`。
2. **scCRAFT（当前安装版）**：`SCIntegrationModel.__init__` 中 `self.device = 'cpu'` 硬编码（CUDA 行被注释），`train_integration_model()` 不接收 device 参数。**所有平台恒走 CPU**，外界无法控制。
3. **Mac MPS 对 scVI/scANVI 的数值稳定性**：未经充分验证，默认保守走 CPU；可通过显式参数覆盖。

如果 notebook 中每个训练 cell 各自判断 `torch.cuda.is_available()` / `torch.backends.mps`，代码分散且容易写错 accelerator 值（scvi-tools 用 `"gpu"` 而非 `"cuda"`，这是常见踩坑点）。

## Decision

### 1. 单点收口：`platform.detect_device()`

在 `src/scrna_integration/platform.py` 新增 `detect_device(prefer, for_method)` 函数，作为框架唯一设备检测入口。遵循 ADR-0010 的 OS 单点收口原则，将设备检测也收口到 `platform.py`。

返回 dict（直接可喂 scvi train）：
```python
{
    "accelerator": str,   # "gpu" | "cpu" | "mps"
    "devices": object,    # "auto" | 1 | [0]
    "device_str": str,    # "cuda" | "mps" | "cpu" (torch 风格，日志用)
    "reason": str,        # 中文决策说明
}
```

决策链：`try import torch` -> CUDA 可用 > MPS 可用 > CPU。`prefer` 参数覆盖自动检测；`for_method` 允许方法级默认行为（scVI/scANVI 在 Mac 上默认 CPU）。

### 2. env_check 集成设备诊断

`platform.env_check()` 在 torch / scvi-tools 版本检测后追加设备检测打印，return dict 新增 `"device"` 字段（`str | None`）。

### 3. 04_embedded notebook 五处改造

| 位置 | 改动 |
|------|------|
| PARAMS cell | 新增 `DEVICE = "auto"` 参数 |
| setup cell | `env_check()` 后打印设备概览 |
| scVI train cell | `model.train()` 加 `accelerator` / `devices` |
| scANVI cell | 两处 `train()` 加 `accelerator` / `devices` |
| scCRAFT cell | 注释说明硬编码 CPU + 日志打印（不传参） |

### 4. 不动的部分

- **`environment.yml` torch pin 不碰**：`torch==2.12.0`（无后缀），pip 在各平台自动装对应 variant（Mac 装 `torch-2.12.0-cp311-none-macosx_11_0_arm64.whl`，Linux 装 `torch-2.12.0-cp311-cp311-manylinux1_x86_64.whl` 等）。
- **其他 src 模块、其他 notebook 不改**。

## 实测约束（后人最易踩，务必读）

1. **scvi-tools 1.4.2** 用 `accelerator`（`"gpu"` / `"cpu"` / `"mps"`）+ `devices`（int / list / `"auto"`），**不是 `use_gpu`**。CUDA 对应 `"gpu"` 不是 `"cuda"`。
2. **scVI / scANVI 在 Mac 上默认 CPU**：`detect_device("auto", for_method="scvi")` 在 MPS 可用时返回 `"cpu"`，因为 MPS 数值稳定性未经充分验证。如需 MPS 实验，显式设 `DEVICE="mps"`。
3. **scCRAFT 恒 CPU**：当前安装版源码 `SCIntegrationModel.__init__` 中 `self.device = 'cpu'` 硬编码。`train_integration_model()` API 不接收 device 参数。`detect_device(for_method="sccraft")` 明确报告此约束。如将来 scCRAFT 上游支持 device 参数，更新本函数的分支逻辑即可。
4. **torch pin**：`environment.yml` 中 `torch==2.12.0`（无 build suffix），pip 自动解析不同平台的 wheel variant。不要手动加 CUDA/MPS suffix——那会破坏跨平台兼容性。

## Considered Options

- **在 notebook 中各自判断 device**：拒绝。分散的 `torch.cuda.is_available()` 调用容易写错 accelerator 值，且背离 ADR-0010 的单点收口原则。
- **环境变量控制**（`DEVICE=cuda`）：辅助方案——`DEVICE` 作为 notebook PARAMS 的默认值存在，`detect_device(prefer=DEVICE)` 读取。PI 可在 notebook 上方改一行来覆盖。
- **完全统一到 CPU**：拒绝。Ubuntu-CUDA 生产环境的存在本身就是为了加速训练，统一 CPU 等于放弃了这台机器的价值。

## Consequences

- `src/scrna_integration/platform.py` 新增 `detect_device()`（~100 行）；`env_check()` return dict 新增 `"device"` key。
- `tests/test_platform.py` 新增 8 个 `detect_device` 测试（mock torch CUDA / MPS / ImportError 分支）。
- `notebooks/04_embedded.ipynb` 5 个 cell 改造：PARAMS 加 `DEVICE`、setup 加概览打印、scVI/scANVI train 加 accelerator/devices、scCRAFT 加注释。
- ADR-0010 "GPU 分歧不存在" 句后追加注记，指向本 ADR。
- code-reviewer 红线更新：notebook/src 中出现 `torch.cuda.is_available()` / `torch.backends.mps` 直接判断 → flag（应通过 `detect_device()` 收口）。
