#!/bin/bash
# Replace CPU torch with CUDA variant on Linux GPU machines.
#
# Run ONCE after: conda env create -f environment.yml && conda activate scrna-integration
# NOT needed on Mac — MPS is built into the standard torch wheel.
# See ADR-0010 + ADR-0013.

set -euo pipefail

ENV_NAME="${1:-scrna-integration}"
CUDA_INDEX="https://download.pytorch.org/whl/cu126"

echo "=== Reinstalling PyTorch 2.12.0 with CUDA 12.6 support ==="
echo "Environment: ${ENV_NAME}"

conda run -n "${ENV_NAME}" pip install torch==2.12.0 \
    --index-url "${CUDA_INDEX}" \
    --force-reinstall \
    --no-deps
# --no-deps 跳过 nvidia-cublas-cu12 等 CUDA 运行时依赖；多数 GPU 集群
# 通过系统 CUDA toolkit 提供（LD_LIBRARY_PATH）。若验证 CUDA not available，
# 手动补装：pip install --index-url "${CUDA_INDEX}" nvidia-cublas-cu12 nvidia-cudnn-cu12

echo ""
echo "=== Verifying CUDA ==="
conda run -n "${ENV_NAME}" python -c "
import torch
print(f'PyTorch: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA version: {torch.version.cuda}')
    print(f'GPU count: {torch.cuda.device_count()}')
    for i in range(torch.cuda.device_count()):
        print(f'  GPU {i}: {torch.cuda.get_device_name(i)}')
else:
    print('WARNING: CUDA not available.')
"
echo "=== Done ==="
