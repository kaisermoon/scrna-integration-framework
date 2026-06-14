# Mac 环境同步清单

> 最近一次 Linux 环境调整：2026-06-13（见 ADR-0012）
> 本文档记录在 Linux 上做过、但 Mac 上尚未同步的环境变更。
> 下次在 Mac 运行前，按此清单同步。

## 背景

为集成 scCRAFT（04_embedded 的可选嵌入方法），Linux 主环境做了调整。Mac 需要同步以下变更才能跑 scCRAFT，并避免 scVI/scCRAFT 训练 segfault。

## Mac 上需要执行的同步

### 1. 主环境 scrna-integration：移除 TensorFlow（必须）

TF 的 oneDNN/MKL 与 PyTorch backward 在同进程冲突会导致训练 segfault。

```bash
conda activate scrna-integration
pip uninstall -y tensorflow tensorflow-probability keras
```

验证：
```bash
python -c "from scrna_integration.platform import env_check; env_check()"
# 应显示 ✓ 无 tensorflow 警告
```

### 2. 安装 scCRAFT（如需用 scCRAFT 嵌入）

```bash
conda activate scrna-integration
cd /tmp
git clone https://github.com/ch2343/scCRAFT
cd scCRAFT && pip install .
# 装完确认 TF 仍未被 scCRAFT 依赖拉回：
pip uninstall -y tensorflow tensorflow-probability keras 2>/dev/null || true
```

### Mac 上的设备支持（已明确）

- **scVI / scANVI**：在 Mac 上**默认走 CPU**。MPS 加速的数值稳定性未经充分验证，框架选择保守策略。如需尝试 MPS，在 `04_embedded` 的 `PARAMS` 中设 `DEVICE="mps"` 即可覆盖。
- **scCRAFT**：当前安装版 `SCIntegrationModel.__init__` 中 `self.device = 'cpu'` 硬编码（CUDA 行被注释），`train_integration_model()` 不接收 device 参数。**所有平台恒走 CPU，Mac MPS 不可达**。
- 详见 ADR-0013。

### 3. scCODA 独立环境（如需用 11_abundance 的 scCODA）

Linux 上已建 `scrna-sccoda` 环境。Mac 上需重新创建（`environment-sccoda.yml` 当前是 linux-64 专属，pin 了 linux 底层库）：

```bash
# Mac 上不能直接 conda env create -f environment-sccoda.yml（含 linux-64 专属库）
# 改为手动创建：
conda create -n scrna-sccoda python=3.11 -y
conda activate scrna-sccoda
pip install sccoda
```

### 4. 刷新环境快照

同步完成后，在 Mac 上重新生成快照，供跨机对比：
```bash
conda activate scrna-integration
python scripts/env_parity.py snapshot   # 生成/更新 docs/env-snapshots/osx-arm64.json
python scripts/env_parity.py compare    # 对比两机差异
```

## 自动检测

每次运行 03-06 主流水线 notebook 时，setup cell 会自动调用 `env_check()`，
检测当前平台 + conda 环境 + TF 是否误入主环境，并打印需要的调整命令。
如看到警告或错误提示，按提示处理。

## 验证同步成功

```bash
conda activate scrna-integration
python -c "
from scrna_integration.platform import env_check
r = env_check()
assert r['ok'], '环境仍有问题'
print('Mac 环境同步完成')
"
```
