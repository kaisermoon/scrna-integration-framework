#!/usr/bin/env bash
# 预拉取 gitleaks pre-commit hook，避免首次 git commit 卡 5-10 秒。
# 使用：项目 kickoff 阶段、或 CI 准备阶段执行一次即可。
# 幂等：已拉取过会直接 skip。

set -euo pipefail

REPO_ROOT="${1:-$(pwd)}"

if [ ! -f "$REPO_ROOT/.pre-commit-config.yaml" ]; then
  echo "[prefetch-gitleaks] no .pre-commit-config.yaml at $REPO_ROOT, skip" >&2
  exit 0
fi

if ! command -v pre-commit >/dev/null 2>&1; then
  echo "[prefetch-gitleaks] pre-commit not installed; run: pip install pre-commit" >&2
  exit 1
fi

# 把 .pre-commit-config.yaml 中所有 hook 的 environment 全部预安装
cd "$REPO_ROOT"
pre-commit install-hooks

echo "[prefetch-gitleaks] OK: hooks env prefetched at $REPO_ROOT" >&2
