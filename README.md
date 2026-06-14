# scRNA-seq整合分析框架

> 跨病种、多源 scRNA-seq 数据整合分析框架。支持 cellranger / h5ad / RData 多源接入，模块化预处理与质控，多方法去批次与注释，下游可扩展（拟时序、GRN）。每步骤可替换、可迭代回跑，以 Jupyter notebook 为主要交付形式。

[![CI](https://github.com/kaisermoon/scrna-integration-framework/actions/workflows/test.yml/badge.svg)](https://github.com/kaisermoon/scrna-integration-framework/actions/workflows/test.yml)
[![Lint](https://github.com/kaisermoon/scrna-integration-framework/actions/workflows/lint.yml/badge.svg)](https://github.com/kaisermoon/scrna-integration-framework/actions/workflows/lint.yml)

## 设计原则

1. **多源接入**：cellranger filtered/raw matrix、h5ad、RData → 统一为 CellxGene-compatible AnnData
2. **模块化与可替换**：每步骤约束输入输出（h5ad），任意环节可换方法做对比
3. **迭代回跑**：注释步骤发现问题 → 回退到 HVG / embedding / clustering 重跑
4. **多方法并存**：去批次（Harmony/scVI/scANVI）、注释（marker/LLM/scANVI/AUCell/UCell）等关键决策步骤允许多方法同时跑，结果存于 obs 列做交叉比对
5. **内存自律**：稀疏矩阵 + 避免对象复制，面向海量整合
6. **学术追新**：新方法发表 → 快速插入流水线，不重写下游

## 仓库结构

```
.
├── src/                  # 模块化 Python 包
│   ├── io/              # 多源数据读取与 schema 标准化
│   ├── qc/              # 质控与过滤
│   ├── preprocessing/   # normalize / HVG
│   ├── embedding/       # PCA / Harmony / scVI / scANVI
│   ├── clustering/      # Leiden / 参数扫描
│   ├── annotation/      # marker / LLM / scANVI / AUCell / UCell
│   └── downstream/      # DEG / 拟时序 / GRN
├── notebooks/           # 编号流水线 notebook（01_io → 09_downstream）
├── scripts/             # 一次性工具脚本
├── tests/               # pytest 单元测试
├── data/                # 外部原始数据（gitignore）
├── results/             # 分析产出（gitignore）
├── references/          # 早期框架蓝本、学生参考代码、标记物库
│   ├── markers/         # 标记物库（PI 维护的基因-细胞类型映射表）
│   │   ├── README.md    #   字段说明与填写指南
│   │   └── *.csv        #   按组织命名的标记物文件（空模板 → PI 填写）
│   ├── legacy-GCPL/     #   PI 历史项目代码
│   └── student-code/    #   学生下游分析技术参考
```

## 依赖

- Python ≥ 3.10
- 核心：`scanpy`、`anndata`、`scvi-tools`、`scvelo`（部分模块）
- R ≥ 4.3（可选）：Monocle3 / InferCNV
- 详见 `pyproject.toml`

## 环境搭建

本项目使用 conda 隔离环境，**禁止装入 base 或主环境**。

```bash
# Python 环境（固定命名 scrna-integration）
conda env create -f environment.yml
conda activate scrna-integration

# scCODA 丰度差异分析环境（固定命名 scrna-sccoda；可选，仅需 scCODA 时安装）
conda env create -f environment-sccoda.yml
conda activate scrna-sccoda

# R 环境（固定命名 scrna-integration-r；可选，仅需下游 R 工具时安装）
conda env create -f environment-r.yml
conda activate scrna-integration-r
```

开发安装（轻量，仅 pytest/ruff/pre-commit）：`conda activate scrna-integration && pip install -e ".[dev]" && pre-commit install`。完整重装流程先确保已激活 `scrna-integration` 环境。

> **环境分离说明**：`scrna-integration`（主环境）为保持 scCRAFT 训练稳定已移除 TensorFlow。scCODA 依赖 TF，因此独立为 `scrna-sccoda`——TF 的 oneDNN/MKL 与 PyTorch 同进程 backward 时存在冲突 segfault。完整分析流程中，scVI/scCRAFT 步骤用主环境，scCODA 丰度分析步骤切换到 `scrna-sccoda`。
>
> **R 重型包单独安装**：Monocle3、hdWGCNA 等无 conda 预编译包，需在激活 `scrna-integration-r` 环境后进入 R 单独安装——`BiocManager::install("monocle3")` / `remotes::install_github("smorabit/hdWGCNA")`。详见 `environment-r.yml` 注释。

## 运行

框架尚处 planning 阶段，暂无 CLI 入口。首期端到端验证 notebook 计划放在 `notebooks/end-to-end-pilot.ipynb`。

## 多环境运行 / 计算设备

项目支持三种环境的计算设备自适应：**Mac (Apple Silicon/MPS)**、**Linux 无显卡 (CPU)**、**Ubuntu (CUDA GPU)**。设备检测收口到 `platform.detect_device()`（见 ADR-0013）。

### 使用方式

`04_embedded` 的 `PARAMS` 中新增 `DEVICE = "auto"`。`auto` 模式自动选择最优设备：CUDA GPU > Mac 下 scVI/scANVI 走 CPU > CPU。也可显式覆盖：

```python
PARAMS = {
    "DEVICE": "auto",     # "auto" | "cuda" | "mps" | "cpu"
    ...
}
```

### 三平台行为

| 环境 | auto 检测结果 | scVI / scANVI | scCRAFT |
|---|---|---|---|
| Ubuntu + CUDA | gpu | GPU 加速 | CPU（内部硬编码） |
| Mac + MPS | cpu（scVI/scANVI） | CPU（MPS 数值稳定性未验证，可显式 `DEVICE='mps'`） | CPU |
| Linux 无显卡 | cpu | CPU | CPU |

> **说明**：scvi-tools 1.4.2 使用 `accelerator`（`"gpu"` / `"cpu"` / `"mps"`）控制设备，非旧参数 `use_gpu`。CUDA 对应 `"gpu"`。scCRAFT 当前安装版 `self.device = 'cpu'` 硬编码，所有平台恒走 CPU。Mac 下如需尝试 MPS，在 PARAMS 设 `DEVICE="mps"` 即可覆盖默认 CPU 策略。

### 环境安装

`environment.yml` 的 torch pin **不带平台后缀**（如 `pytorch=2.12.0=*_0`），各平台 `pip install` 自动安装对应 variant（CUDA / MPS / CPU）。生产部署（Ubuntu CUDA）通过 `git pull` / PR 同步代码，跨机环境对齐用 `scripts/env_parity.py`（见 ADR-0010）。详见 [ADR-0013](docs/adr/0013-device-adaptive-layer.md)。

## 标记物库

框架提供标记物库（`references/markers/`），用于管理细胞类型与标记基因的对应关系。标记物库由 PI（消化科专家）维护，框架**不预置任何真实标记物基因**——所有标记物内容由 PI 亲自填写。

**快速开始**：

```python
from scrna_integration import load_markers

# 加载标记物（默认返回 canonical + optional 标记物）
markers = load_markers("references/markers/gastric_epithelial.csv")

# 使用前检查基因存在性（不同数据集基因集不同——必须检查）
for cell_type, genes in markers.items():
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"[{cell_type}] {len(missing)} genes missing: {missing}")
    markers[cell_type] = present

# 在 dotplot 中使用
sc.pl.dotplot(adata, var_names=markers, groupby="leiden")
```

详细说明见 [`references/markers/README.md`](references/markers/README.md)（字段定义、role 取值规则、PMID 硬要求、面向消化科专家的填写指南）。

## 测试

```bash
pytest tests/ -v
ruff check src/ tests/
```

## 架构决策

框架的重大架构决策以 ADR（Architecture Decision Record）形式记录在 [`docs/adr/`](docs/adr/) 目录中，涵盖为何选择"薄框架"（不包装 scanpy）、为何只有 3 个公开函数、R 互操作如何分流等关键决策。详见 [ADR 索引](docs/adr/index.md)。

## 贡献

PR 流程见 `.github/PULL_REQUEST_TEMPLATE.md`。提交前 `pre-commit` hook 会自动跑 gitleaks 与基础检查。

## License

MIT（待确认）

---

**作者**：钟子劭（kaisermoon）
**关联科研方向**：AI × 中西医结合 / 消化病单细胞研究
**论文引用**：待项目成熟后补 BibTeX 与 software paper 链接
