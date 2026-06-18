"""IO 模块：多来源单细胞数据读取，manifest 驱动的 obs schema 标准化。

唯一公开函数 ``read_with_manifest``，按照 SPEC 的 13 个行为实现。
主体按步骤分块，每块中文注释解释在做什么、为什么。
面向非计算机专业 PI/学生：从上到下线性可读，打开文件就能看懂整个数据读取流程。

复杂逻辑保留为私有 helper（含充分中文注释）：
- _read_10x_mtx / _find_10x_dirs（多目录发现）
- _read_txt_gz（Yue organoid 格式）
- _read_rds（stub，待 R 环境）
- _sync_gene_ids / _sync_symbol_to_ensembl / _mygene_symbol_to_ensembl / _sync_ensembl_to_symbol（基因 ID 双向同步）
- _join_one_clinical_table（临床表关联）
- _warn_layer2（Layer 2 字段校验）
- _compute_baseline_qc（基线 QC 指标）
- _load_manifest / _validate_manifest（manifest 加载与校验）
"""

from __future__ import annotations

import os
from collections import Counter
import warnings
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

# =============================================================================
# 公开入口
# =============================================================================


def read_with_manifest(manifest_path: str) -> anndata.AnnData:
    """读入并标准化单细胞数据——框架核心入口函数。

    按照 manifest YAML 配置完成 13 个标准化步骤（详见 SPEC "read_with_manifest"）。
    读取完成后返回普通 AnnData，调用方可自由写入 h5ad 或继续下游分析。

    面向非计算机专业学生：函数主体按步骤分块，每块 # === 注释 === 解释
    在做什么、为什么这个步骤是必要的。
    """
    # ===== 0. 加载与校验 manifest =====
    manifest = _load_manifest(manifest_path)
    _validate_manifest(manifest)

    source_dataset = str(manifest["source_dataset"])
    input_block = manifest["input"]
    fmt = input_block["format"]
    path = input_block["path"]

    # ===== 1. 读矩阵（按 format 分发）=====
    # 不同数据集的原始存储格式不同：
    #   10x_mtx → cellranger 输出目录（最常见）
    #   h5ad → 别人已处理好的 AnnData 文件
    #   h5 → 10x 官方的 HDF5 格式
    #   txt.gz → 某些研究者自制的压缩文本矩阵
    #   rds → Seurat R 对象（待支持）
    # 为什么不在 manifest 里统一成一种格式：真实科研数据来源多样，
    # 强制转换格式反而增加上游负担；框架做分发比上游做转换更省力。
    if fmt == "h5ad":
        # h5ad 是最简单的：直接 scanpy 读入，补上数据集标签
        adata = sc.read_h5ad(path)
        adata.obs["source_dataset"] = source_dataset
    elif fmt == "h5":
        adata = sc.read_10x_h5(path)
        adata.obs["source_dataset"] = source_dataset
    elif fmt == "10x_mtx":
        # 10x mtx 需要处理多样本目录结构（见 _read_10x_mtx）
        adata = _read_10x_mtx(path, source_dataset)
    elif fmt == "txt.gz":
        # Yue 等非标格式（见 _read_txt_gz）
        adata = _read_txt_gz(path, source_dataset)
    elif fmt == "rds":
        raise NotImplementedError(
            "RDS 格式尚未支持。Tsubosaka 数据集需要 R 环境（ADR-0007）。"
            "详见 docs/adr/0007-r-bridge-tool-split.md。"
        )
    else:
        raise ValueError(
            f"不支持的 input.format '{fmt}'。支持: 10x_mtx, h5ad, h5, txt.gz, rds"
        )

    # ===== 2. obs 列映射 + 值映射 =====
    # 为什么需要这两步：不同数据集的 obs 列名各异——
    # 有的叫 Sample/Patient/Group，有的叫 orig.ident/group/condition。
    # manifest 的 obs_mapping 定义"源列→目标标准列"的对应关系。
    # value_mapping 进一步统一值（如把 CAG_mild / CAG_severe 归为 CAG），
    # 避免下游分析因命名细微差异而分裂分组。
    obs_mapping = manifest.get("obs_mapping", {})
    value_mapping = manifest.get("value_mapping", {})
    if obs_mapping:
        for target_col, source_col in obs_mapping.items():
            if source_col not in adata.obs.columns:
                warnings.warn(
                    f"obs_mapping: 源列 '{source_col}' 在数据中未找到；无法填充 '{target_col}'",
                    stacklevel=2,
                )
                continue

            # P0-2: 目标列冲突检查——已有列被覆盖时静默是隐患
            if target_col in adata.obs.columns:
                warnings.warn(
                    f"obs_mapping: 目标列 '{target_col}' 已存在，将被覆盖",
                    stacklevel=2,
                )

            # DG-4: categorical 列 .astype(str) 会丢失类别排序信息
            if isinstance(adata.obs[source_col].dtype, pd.CategoricalDtype):
                warnings.warn(
                    f"obs_mapping: 源列 '{source_col}' 是 categorical 类型，"
                    f"转换为字符串将丢失类别排序信息",
                    stacklevel=2,
                )

            # P0-4: NaN 保留为真 NaN（而非字符串 "nan"）
            source_series = adata.obs[source_col]
            nan_mask = source_series.isna()
            adata.obs[target_col] = source_series.astype(str)
            if nan_mask.any():
                adata.obs.loc[nan_mask, target_col] = np.nan

        for target_col, mapping in value_mapping.items():
            if target_col not in adata.obs.columns:
                warnings.warn(
                    f"value_mapping: 目标列 '{target_col}' 未找到；跳过",
                    stacklevel=2,
                )
                continue

            # P0-3: 检测未被 mapping 覆盖的取值，不静默穿透
            col_vals = adata.obs[target_col]
            non_null_vals = col_vals.dropna().astype(str)
            unique_vals = non_null_vals.unique()
            uncovered = sorted([v for v in unique_vals if v not in mapping])
            if uncovered:
                # 附带频数，让 PI 可见各源词表分歧
                val_counts = Counter(non_null_vals)
                uncovered_detail = {
                    v: val_counts[v] for v in uncovered
                }
                warnings.warn(
                    f"value_mapping: 目标列 '{target_col}' 有 {len(uncovered)} 个取值"
                    f"未被映射覆盖，将保持原值穿透: {uncovered_detail}",
                    stacklevel=2,
                )

            adata.obs[target_col] = adata.obs[target_col].map(
                lambda x, m=mapping: m.get(str(x), x)
            )

    # ===== 3. 关联临床信息表 =====
    # 外部临床 metadata（xlsx/csv）中的 patient/sample 级字段
    # （如年龄、性别、H. pylori 状态）需要合并到 obs 中。
    # 每个表通过 join_on 指定关联键，on_missing/on_conflict 控制冲突策略。
    tables = manifest.get("clinical_metadata", [])
    if tables:
        for tbl_cfg in tables:
            _join_one_clinical_table(adata, tbl_cfg)

    # ===== 4. 注入本体论与项目级常量 =====
    # ontology: 该数据集固定的本体信息（如组织="gastric mucosa"），
    #   每行细胞都写入相同的值，确保跨数据集分析时这些维度统一可用。
    # project_specific: 项目自定义字段，可基于已有 obs 列映射
    #   （如根据 Group 列推导 disease_grade），也可直接赋常量值。
    ontology = manifest.get("ontology", {})
    for key, value in ontology.items():
        adata.obs[key] = value

    project_specific = manifest.get("project_specific", {})
    for col_name, cfg in project_specific.items():
        if isinstance(cfg, dict) and "source_column" in cfg:
            src = cfg["source_column"]
            rules = cfg.get("rules", {})
            if src in adata.obs.columns:
                adata.obs[col_name] = adata.obs[src].map(
                    lambda x, r=rules: r.get(str(x), str(x))
                )
            else:
                # P0-6: 源列缺失时明确警告，消除全文件唯一静默跳过路径
                warnings.warn(
                    f"project_specific: 源列 '{src}' 在数据中未找到；"
                    f"无法填充 '{col_name}'",
                    stacklevel=2,
                )
        elif isinstance(cfg, dict) and "value" in cfg:
            adata.obs[col_name] = cfg["value"]

    # ===== 5. 生成全局唯一细胞 ID =====
    # 为什么要自己生成：多个数据集整合后 barcode 会冲突。
    # 例如 AAACCCAAGAA-1 可能同时存在于 GSE249874 和 GSE134520 中，
    # 但它们是两个完全不同的细胞。cell_id = {数据集}_{样本}_{barcode} 确保唯一。
    # barcode 去掉 -1 后缀（那是 10x GEM well suffix，不是细胞身份的一部分）。
    sample = adata.obs.get("sample_id", pd.Series("unknown", index=adata.obs_names))
    barcode = [str(bc).split("-")[0] for bc in adata.obs_names]
    adata.obs["cell_id"] = [
        f"{source_dataset}_{s}_{b}" for s, b in zip(sample, barcode, strict=True)
    ]

    # ===== 6. 原始作者标注列重命名 =====
    # 格式: cell_type_original_{数据集名}_v1[{_角色}]
    # 为什么需要重命名：不同数据集的原始作者标注列名各异
    # （cell_type / CellType / cluster_name...），合并后需要统一命名模式
    # 以便 stage 6 交叉比对时按来源追踪标注差异。
    annotations = manifest.get("original_annotations", [])
    for entry in annotations:
        col = entry["column"]
        role = entry.get("role", "")
        if col not in adata.obs.columns:
            warnings.warn(
                f"original_annotations: 列 '{col}' 在数据中未找到；跳过",
                stacklevel=2,
            )
            continue
        suffix = f"_{role}" if role else ""
        new_name = f"cell_type_original_{source_dataset}_v1{suffix}"
        adata.obs[new_name] = adata.obs[col]

    # ===== 7. 基因 ID 双向同步 =====
    # 确保 var.index = gene symbol（符合 scanpy/scverse 生态惯例），
    # 同时 var["ensembl_id"] 存储 Ensembl ID（供跨数据库查询、mygene 转换等）。
    # 为什么需要双向：scanpy 的所有函数默认用 var.index 做基因名查找，
    # 但跨物种比较、在线查询、某些 atlas 对齐需要 Ensembl ID。
    # 为什么不用 Ensembl 做主索引：学生/合作者阅读 h5ad 时 human-readable
    # 的 gene symbol 远比 ENSG 编号友好，且 scanpy 生态默认符号。
    _sync_gene_ids(adata, input_block.get("gene_id_format", "auto"))

    # ===== 8. 物种校验 =====
    # 目前仅接受 human。为什么：不同物种的基因命名体系（小鼠首字母大写）、
    # 标记物参考文献、参考基因组都完全不同。支持跨物种需要重新设计
    # markers 库、QC 阈值、mygene 查询逻辑——那是另一条技术路线。
    species = manifest["species"]
    if species != "human":
        raise ValueError(
            f"物种 '{species}' 不受支持。框架目前仅接受 species='human'。"
            f"其他物种需要新 ADR 论证技术路线。"
        )
    adata.uns["species"] = "human"

    # ===== 9. 疾病系统 + 项目 ID 传播 =====
    # disease_system / project_id 是跨项目细胞分组的 Layer 1 必需字段。
    # 每行细胞标记所属系统与项目，使跨疾病整合分析能按系统分组。
    adata.obs["disease_system"] = manifest["disease_system"]
    adata.obs["project_id"] = manifest["project_id"]

    # ===== 10. Layer 2 强警告 =====
    # CellxGene 对齐的 7 个字段缺失或异常时发出警告但不阻断读取。
    # LLM best-effort 修复延后到 PR-3c（需要 OpenRouter key）。
    _warn_layer2(adata)

    # P0-1: Layer 1 确定性校验——三个必需字段缺任一即 fail loudly
    _validate_layer1(adata)

    # ===== 11. 记录 raw matrix 路径 =====
    # 为什么只记路径不加载：raw matrix（未过滤的 cellranger 输出）通常是
    # filtered matrix 的 2-3 倍大。stage 1 加载会造成内存翻倍，而只有
    # stage 2 的 SoupX 环境校正才真正需要它。按需加载更安全。
    raw_path = input_block.get("raw_path")
    adata.uns["raw_matrix_path"] = raw_path if raw_path else None

    # ===== 12. 计算基线 QC 指标 =====
    # 直接在 adata.X（原始 counts）上计算 n_genes / total_counts /
    # pct_counts_mt / pct_counts_ribo。为什么放在读取阶段：
    # 确保无论上游数据是否已做预处理，这些列在每个数据集中都已存在且对齐，
    # 避免 stage 2 QC 遇到"这个数据集的 pct_mt 列去哪了"的问题。
    _compute_baseline_qc(adata)

    # ===== 13. 返回 AnnData =====
    # 返回的是普通 AnnData 对象，没有框架特有的元数据封装。
    # 调用方可以自由地 write_h5ad 做 checkpoint，或继续 scanpy 下游操作。
    return adata


# =============================================================================
# 10x mtx 读取（保留，复杂：多样本目录发现 + 拼接）
# =============================================================================


def _read_10x_mtx(
    path: str, source_dataset: str, _input_block: dict | None = None
) -> anndata.AnnData:
    """读取一个或多个 10x mtx 目录（cellranger filtered_feature_bc_matrix）。

    *path* 可以是单个 ``filtered_feature_bc_matrix`` 目录，
    也可以是包含多个子样本目录的父目录（如 Nancang 数据集的批量结构）。

    为什么不用 scanpy 内置的多样本读取：scanpy 没有"在父目录下自动发现
    所有 filtered_feature_bc_matrix 子目录并拼接"的功能——每个样本跑一次
    sc.read_10x_mtx 然后 concat。这个 helper 封装了这段机械化操作。
    """
    mtx_path = Path(path)
    subdirs = _find_10x_dirs(mtx_path)
    if not subdirs:
        raise FileNotFoundError(
            f"在 {path} 下未找到 10x mtx 目录"
        )

    parts: list[anndata.AnnData] = []
    for sub in sorted(subdirs):
        # 样本名：如果子目录叫 filtered_feature_bc_matrix，取父目录名
        sample_name = (
            sub.parent.name if sub.name == "filtered_feature_bc_matrix" else sub.name
        )
        try:
            adata_part = sc.read_10x_mtx(sub, var_names="gene_symbols", cache=False)
        except Exception:
            # 某些 10x 目录可能没有 gene_symbols；回退到 gene_ids
            adata_part = sc.read_10x_mtx(sub, var_names="gene_ids", cache=False)
        adata_part.obs["source_dataset"] = source_dataset
        adata_part.obs["sample_id"] = sample_name
        # barcode 加样本前缀防冲突
        adata_part.obs_names = [f"{sample_name}_{bc}" for bc in adata_part.obs_names]
        parts.append(adata_part)

    if len(parts) == 1:
        adata = parts[0]
    else:
        adata = anndata.concat(parts, join="outer", index_unique="-")
        # sc.read_10x_mtx 对缺失基因自动填 0，但 concat(join='outer') 可能
        # 产生 NaN 并 densify 矩阵；强制转回 sparse 并填 0
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(np.nan_to_num(adata.X, nan=0))

    return adata


def _find_10x_dirs(root: Path) -> list[Path]:
    """在 *root* 下找到所有 filtered_feature_bc_matrix 目录。

    三层搜索策略：
    1. root 本身就是 mtx 目录（含有 matrix.mtx.gz 或 matrix.mtx）
    2. root 下直接有 filtered_feature_bc_matrix 子目录
    3. root 的每个子目录下找 filtered_feature_bc_matrix（多样本批量结构）
    """
    # root 本身就是目标
    if (root / "matrix.mtx.gz").exists() or (root / "matrix.mtx").exists():
        return [root]
    # root 下有 filtered_feature_bc_matrix
    if (root / "filtered_feature_bc_matrix").is_dir():
        return [root / "filtered_feature_bc_matrix"]
    # 搜索子目录
    result = []
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        candidate = child / "filtered_feature_bc_matrix"
        if candidate.is_dir():
            result.append(candidate)
    return result


# =============================================================================
# txt.gz 读取（保留，复杂：非标格式，基因×细胞需转置）
# =============================================================================


def _read_txt_gz(
    path: str, source_dataset: str, _input_block: dict | None = None
) -> anndata.AnnData:
    """读取制表符分隔的 ``*.txt.gz`` 计数文件（Yue organoid 格式）。

    每个文件：第一列 = 基因，其余列 = 细胞 barcode。
    矩阵是基因×细胞，加载后转置为细胞×基因（scanpy 标准方向）。
    为什么保留这个 reader：这种格式在非标准数据集中常见，
    且转置 + sparse 转换 + 多文件拼接的逻辑不琐碎。
    """
    txt_dir = Path(path)
    files = sorted(txt_dir.glob("*.txt.gz"))
    if not files:
        raise FileNotFoundError(f"在 {path} 下未找到 *.txt.gz 文件")

    parts: list[anndata.AnnData] = []
    for fp in files:
        sample_id = fp.stem.replace("_count", "").split(".")[0]
        # 读入 tab 分隔的稠密矩阵（基因×细胞）
        df = pd.read_csv(fp, sep="\t", index_col=0, compression="gzip")
        # 转置为细胞×基因（scanpy 标准方向）
        X = sp.csr_matrix(df.values.T.astype(np.float32))  # noqa: N806
        var = pd.DataFrame(index=df.index.values)
        obs = pd.DataFrame(index=df.columns.values)
        obs["source_dataset"] = source_dataset
        obs["sample_id"] = sample_id
        adata_part = anndata.AnnData(X=X, obs=obs, var=var)
        parts.append(adata_part)

    if len(parts) == 1:
        adata = parts[0]
    else:
        adata = anndata.concat(parts, join="outer", index_unique="-")
        if not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(np.nan_to_num(adata.X, nan=0))

    return adata


def _read_rds(
    _path: str, _source_dataset: str, _input_block: dict | None = None
) -> anndata.AnnData:
    """RDS 读取器——**尚未实现**（需要 R 环境，见 ADR-0007）。

    Tsubosaka 数据集是 RDS 格式，在 PR-1 夹具中暂被跳过。
    后续 PR 将通过 rpy2 实现 RDS 读取。
    """
    raise NotImplementedError(
        "RDS 格式尚未支持。"
        "Tsubosaka 数据集需要 R 环境（ADR-0007）。"
        "详见 docs/adr/0007-r-bridge-tool-split.md。"
        "将在后续 PR 中实现。"
    )


# =============================================================================
# 临床 metadata 关联（保留，复杂：冲突策略、左连接、重命名）
# =============================================================================


def _join_one_clinical_table(adata: anndata.AnnData, cfg: dict) -> None:
    """将一个临床信息表格关联到 adata.obs（左连接）。

    为什么保留：涉及多种文件格式（csv/xlsx）、skip_rows、
    column_mapping 重命名、value_mapping 值转换、on_missing/on_conflict
    多重策略控制——内联到 read_with_manifest 会让主函数膨胀得难以阅读。
    """
    file_path = cfg["file"]
    sheet = cfg.get("sheet", 0)
    skip_rows = cfg.get("skip_rows", 0)
    join_on = cfg["join_on"]
    col_mapping = cfg.get("column_mapping", {})
    val_mapping = cfg.get("value_mapping", {})
    on_missing = cfg.get("on_missing", "warn")
    # on_conflict reserved for future use (SPEC "metadata_wins / obs_wins / error")

    manifest_field = join_on["manifest_field"]
    table_column = join_on["table_column"]

    # 文件不存在时的策略
    if not os.path.exists(file_path):
        if on_missing == "strict":
            raise FileNotFoundError(f"临床 metadata 文件未找到: {file_path}")
        warnings.warn(f"临床 metadata 文件未找到: {file_path}", stacklevel=2)
        return

    # 读入（csv 或 xlsx）
    if file_path.endswith(".csv"):
        meta_df = pd.read_csv(file_path)
    else:
        meta_df = pd.read_excel(file_path, sheet_name=sheet, skiprows=skip_rows)

    # 列名映射：把外部表的列名统一到框架标准名
    if col_mapping:
        meta_df = meta_df.rename(columns={v: k for k, v in col_mapping.items()})

    # 值映射：统一取值（如 M→male, F→female）
    for col, mapping in val_mapping.items():
        if col in meta_df.columns:
            meta_df[col] = meta_df[col].map(lambda x, m=mapping: m.get(str(x), x))

    # 检查关联键是否存在于 obs
    if manifest_field not in adata.obs.columns:
        warnings.warn(
            f"临床关联: manifest 字段 '{manifest_field}' 在 obs 中未找到；"
            f"无法关联 {file_path}",
            stacklevel=2,
        )
        return

    # P0-5: 显式类型对齐——避免 category-vs-int64 静默零匹配 / object-vs-int64 崩溃
    n_cells_before = adata.n_obs
    if manifest_field in adata.obs.columns:
        adata.obs[manifest_field] = adata.obs[manifest_field].astype(str)
    if table_column in meta_df.columns:
        meta_df[table_column] = meta_df[table_column].astype(str)

    # 左连接：obs（左）← metadata（右）
    original_cols = set(adata.obs.columns)
    meta_cols_to_merge = [table_column] + [
        c for c in meta_df.columns
        if c != table_column and c not in original_cols
    ]
    adata.obs = adata.obs.reset_index().merge(
        meta_df[meta_cols_to_merge],
        how="left",
        left_on=manifest_field,
        right_on=table_column,
    ).set_index("index")
    adata.obs.index.name = None

    # 行数不应因左连接变化（细胞数不变定律）
    if adata.n_obs != n_cells_before:
        warnings.warn(
            f"临床关联后细胞数变化: {n_cells_before} → {adata.n_obs}。"
            f"键类型不一致（category/int64/object）可能导致静默零匹配，"
            f"已自动执行 astype(str) 对齐。若仍不一致，请检查 join_on 配置。",
            stacklevel=2,
        )


# =============================================================================
# 基因 ID 双向同步（保留，复杂：mygene 在线查询 + 批量处理）
# =============================================================================


def _sync_gene_ids(adata: anndata.AnnData, gene_id_format: str = "auto") -> None:
    """确保 var.index 是 gene symbol 且 var['ensembl_id'] 有 Ensembl ID。

    判断逻辑：
    - var.index 是 symbol（非 ENSG 开头）→ 尝试从现有列提取或 mygene 查询 Ensembl ID
    - var.index 是 ENSG → 尝试转换为 symbol 并设为 index，原 ENSG 存入 var['ensembl_id']
    - gene_id_format 参数可覆盖自动检测

    为什么需要这一步：scanpy 生态默认 var.index = gene symbol，
    但某些数据集（如 CELLxGENE Census 导出的 h5ad）var.index 是 Ensembl ID，
    不转换会导致 sc.pp.normalize_total 等函数因为 key lookup 失败而 crash。
    """
    idx_sample = str(adata.var.index[0])
    is_ensembl = idx_sample.startswith("ENSG")

    if gene_id_format == "symbol" or (gene_id_format == "auto" and not is_ensembl):
        _sync_symbol_to_ensembl(adata)
    elif gene_id_format == "ensembl" or (gene_id_format == "auto" and is_ensembl):
        _sync_ensembl_to_symbol(adata)
    else:
        _sync_symbol_to_ensembl(adata)


def _sync_symbol_to_ensembl(adata: anndata.AnnData) -> None:
    """var.index 是 symbol；尝试添加 var['ensembl_id']。

    优先使用已存在的列（10x 惯例的 gene_ids 列），不存在时回退 mygene 在线查询。
    """
    # 已有 ensembl_id 列且非全空 → 不用重复查询
    if "ensembl_id" in adata.var.columns and adata.var["ensembl_id"].notna().any():
        return

    # 检查 10x 惯例的 gene_ids 列
    if "gene_ids" in adata.var.columns:
        gene_ids = adata.var["gene_ids"]
        adata.var["ensembl_id"] = gene_ids.where(
            gene_ids.astype(str).str.startswith("ENSG"), ""
        )
        n_mapped = (adata.var["ensembl_id"] != "").sum()
        if n_mapped > 0:
            return

    # 回退：在线查询 mygene
    _mygene_symbol_to_ensembl(adata)


def _mygene_symbol_to_ensembl(adata: anndata.AnnData) -> None:
    """通过 mygene.info 在线查询 gene symbol → Ensembl ID。

    批量查询（每批 1000 个基因），空结果留空并在末尾汇总警告。
    为什么需要在线查询：本地没有完整的 symbol↔Ensembl 映射表，
    mygene 是生物信息学界最常用的基因 ID 转换 API。
    """
    import mygene

    mg = mygene.MyGeneInfo()
    symbols = list(adata.var.index)
    ensembl_ids: dict[str, str] = {}

    batch_size = 1000
    for i in range(0, len(symbols), batch_size):
        batch = symbols[i: i + batch_size]
        try:
            results = mg.querymany(
                batch, scopes="symbol", fields="ensembl.gene",
                species="human", returnall=False,
            )
        except Exception:
            warnings.warn(
                f"mygene 查询失败（批次 {i}）；对应基因的 ensembl_id 留空",
                stacklevel=2,
            )
            continue
        for r in results:
            query = r.get("query", "")
            ensembl_data = r.get("ensembl")
            if ensembl_data and isinstance(ensembl_data, dict):
                eid = ensembl_data.get("gene")
                if eid:
                    ensembl_ids[query] = eid

    adata.var["ensembl_id"] = adata.var.index.map(ensembl_ids).fillna("")
    n_missing = (adata.var["ensembl_id"] == "").sum()
    if n_missing > 0:
        warnings.warn(
            f"基因 ID 同步: {n_missing}/{len(symbols)} 个 symbol 无法通过 "
            f"mygene 映射到 Ensembl ID。对应基因的 ensembl_id 留空。",
            stacklevel=2,
        )


def _sync_ensembl_to_symbol(adata: anndata.AnnData) -> None:
    """var.index 是 Ensembl ID；转换为 gene symbol。

    优先使用 var 中的 feature_name 列（CELLxGENE 惯例），
    不存在则回退 mygene 在线查询。
    """
    # CELLxGENE 导出的 h5ad 常有 feature_name 列存储 gene symbol
    if "feature_name" in adata.var.columns:
        symbol_map = adata.var["feature_name"].to_dict()
        adata.var["_orig_ensembl"] = adata.var.index.values
        new_index = [symbol_map.get(eid, eid) for eid in adata.var.index]
        adata.var.index = new_index
        adata.var.index.name = None
        adata.var["ensembl_id"] = adata.var["_orig_ensembl"]
        del adata.var["_orig_ensembl"]
        n_mapped = sum(1 for v in new_index if not str(v).startswith("ENSG"))
        if n_mapped > 0:
            return
        # feature_name 没帮上忙（可能也是 ENSG 格式）；回退 mygene

    # 在线查询 mygene：Ensembl ID → gene symbol
    import mygene

    mg = mygene.MyGeneInfo()
    ensembl_ids = list(adata.var.index)
    symbol_map_result: dict[str, str] = {}

    batch_size = 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i: i + batch_size]
        try:
            results = mg.querymany(
                batch, scopes="ensembl.gene", fields="symbol",
                species="human", returnall=False,
            )
        except Exception:
            warnings.warn(
                f"mygene 查询失败（批次 {i}）；保留 Ensembl ID 作为索引",
                stacklevel=2,
            )
            continue
        for r in results:
            query = r.get("query", "")
            sym = r.get("symbol")
            if sym:
                symbol_map_result[query] = sym

    adata.var["_orig_ensembl"] = adata.var.index.values
    new_index = [symbol_map_result.get(eid, eid) for eid in adata.var.index]
    adata.var.index = new_index
    adata.var.index.name = None
    adata.var["ensembl_id"] = adata.var["_orig_ensembl"]
    del adata.var["_orig_ensembl"]

    n_missing = sum(1 for v in new_index if str(v).startswith("ENSG"))
    if n_missing > 0:
        warnings.warn(
            f"基因 ID 同步: {n_missing}/{len(ensembl_ids)} 个 Ensembl ID 无法通过 "
            f"mygene 映射到 gene symbol。保留 Ensembl ID 作为索引，"
            f"对应基因的 ensembl_id 留空。",
            stacklevel=2,
        )
        mask = adata.var.index.astype(str).str.startswith("ENSG")
        adata.var.loc[mask, "ensembl_id"] = ""


# =============================================================================
# 基因组位置注入（公开，用于 CNV 推断等下游分析）
# =============================================================================


def inject_genomic_positions(
    adata: anndata.AnnData,
    species: str = "human",
    batch_size: int = 1000,
) -> anndata.AnnData:
    """向 adata.var 注入基因组位置信息（chromosome/start/end）。

    使用 mygene.info 的 genomic_pos 端点批量查询。
    查询失败的基因对应列填 NaN 并打印 warning（不中断）。

    Parameters
    ----------
    adata : AnnData
        输入对象。var 的 index 或 'ensembl_id'/'gene_ids' 列作为查询键。
    species : str
        物种，默认 "human"。
    batch_size : int
        每次 API 请求的基因数，默认 1000。

    Returns
    -------
    AnnData
        同一对象（inplace 修改 var，新增 'chromosome'/'start'/'end' 三列）。
    """
    import re

    try:
        import mygene
    except ImportError:
        raise ImportError(
            "mygene 未安装。请运行: pip install mygene"
        ) from None

    # ---- 自动检测基因标识列 ----
    if "ensembl_id" in adata.var.columns:
        gene_ids_raw = adata.var["ensembl_id"].astype(str).tolist()
    elif "gene_ids" in adata.var.columns:
        gene_ids_raw = adata.var["gene_ids"].astype(str).tolist()
    else:
        gene_ids_raw = adata.var.index.astype(str).tolist()

    # ---- 清洗：去掉 Ensembl 版本号后缀 ----
    gene_ids_clean = [re.sub(r"\.\d+$", "", str(g)) for g in gene_ids_raw]

    # ---- 推断查询 scope ----
    _sample_ids = gene_ids_clean[: min(5, len(gene_ids_clean))]
    if any(str(g).startswith("ENS") for g in _sample_ids):
        scopes = "ensembl.gene"
    else:
        scopes = "symbol"

    # ---- 批量查询 mygene ----
    mg = mygene.MyGeneInfo()
    pos_map: dict[str, dict] = {}  # clean_id -> {chromosome, start, end}
    n_total = len(gene_ids_clean)

    for i in range(0, n_total, batch_size):
        batch = gene_ids_clean[i : i + batch_size]
        try:
            results = mg.querymany(
                batch,
                scopes=scopes,
                fields="genomic_pos",
                species=species,
                returnall=False,
            )
        except Exception:
            warnings.warn(
                f"mygene 查询失败（批次 {i}）；对应基因的位置信息留空",
                stacklevel=2,
            )
            continue

        for r in results:
            query = r.get("query", "")
            if r.get("notfound"):
                continue
            gpos = r.get("genomic_pos")
            if gpos is None:
                continue
            if isinstance(gpos, list):
                if len(gpos) == 0:
                    continue
                gpos = gpos[0]
            if not isinstance(gpos, dict):
                continue

            chr_val = str(gpos.get("chr", ""))
            if not chr_val:
                continue
            # chr 前缀归一化
            chr_val = re.sub(r"^chr", "", chr_val, flags=re.IGNORECASE)

            pos_map[query] = {
                "chromosome": chr_val,
                "start": gpos.get("start"),
                "end": gpos.get("end"),
            }

    # ---- 写入 adata.var ----
    adata.var["chromosome"] = [
        pos_map.get(gid, {}).get("chromosome", np.nan) for gid in gene_ids_clean
    ]
    adata.var["start"] = pd.to_numeric(
        [pos_map.get(gid, {}).get("start", np.nan) for gid in gene_ids_clean],
        errors="coerce",
    )
    adata.var["end"] = pd.to_numeric(
        [pos_map.get(gid, {}).get("end", np.nan) for gid in gene_ids_clean],
        errors="coerce",
    )

    # ---- 汇总警告 ----
    n_failed = adata.var["chromosome"].isna().sum()
    if n_failed > 0:
        warnings.warn(
            f"⚠️ {n_failed}/{n_total} 基因未获取到位置信息",
            stacklevel=2,
        )

    return adata


# =============================================================================
# Layer 2 警告（保留，规则逻辑不琐碎）
# =============================================================================

_LAYER2_FIELDS = [
    "disease",
    "disease_ontology_term_id",
    "tissue",
    "tissue_ontology_term_id",
    "assay",
    "sex",
    "development_stage",
]


def _warn_layer2(adata: anndata.AnnData) -> None:
    """对缺失/异常的 Layer 2 CellxGene 对齐字段发出强警告。

    LLM best-effort 修复延后到 PR-3c（需要 OpenRouter key）。
    本函数仅警告并留 NaN——绝不阻断读取。
    """
    for field in _LAYER2_FIELDS:
        if field not in adata.obs.columns:
            warnings.warn(
                f"[Layer2] obs 列 '{field}' 缺失。"
                f"LLM best-effort 修复延后到 PR-3c。"
                f"所有细胞将写入 NaN。",
                stacklevel=2,
            )
            adata.obs[field] = np.nan
        else:
            col = adata.obs[field]
            n_null = col.isna().sum()
            n_empty = (col.astype(str).str.strip() == "").sum()
            n_problem = n_null + n_empty
            if n_problem > 0:
                warnings.warn(
                    f"[Layer2] obs 列 '{field}' 有 {n_problem}/{len(col)} 个"
                    f"缺失或空值。LLM best-effort 修复延后到 PR-3c。",
                    stacklevel=2,
                )
            # 检查"unknown""NA"等可疑值
            suspicious = col.astype(str).str.lower().isin(
                ["unknown", "na", "n.a.", "n/a", "none", "null"]
            )
            if suspicious.any():
                warnings.warn(
                    f"[Layer2] obs 列 '{field}' 有 {suspicious.sum()} 个细胞的值"
                    f"是可疑占位符（'unknown'/'NA' 等）。"
                    f"LLM 修复延后到 PR-3c。",
                    stacklevel=2,
                )


# =============================================================================
# Layer 1 确定性校验（P0-1：fail loudly）
# =============================================================================

_LAYER1_REQUIRED = ["source_dataset", "project_id", "disease_system"]


def _validate_layer1(adata: anndata.AnnData) -> None:
    """校验 Layer 1 三个必需字段是否已写入 adata.obs。

    缺任一即 raise ValueError（fail loudly，PI 已拍板抛错不是 warn）。
    这是纯确定性校验——零 LLM，只验证 manifest 冻结的映射是否已正确应用。
    """
    missing = [f for f in _LAYER1_REQUIRED if f not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Layer 1 必需字段缺失: {missing}。"
            f"read_with_manifest 未正确写入这些字段，请检查 manifest 配置"
            f"（source_dataset/project_id/disease_system）"
            f"及 obs_mapping 映射是否正确。"
        )


# =============================================================================
# 基线 QC 指标（保留，计算逻辑不琐碎）
# =============================================================================


def _compute_baseline_qc(adata: anndata.AnnData) -> None:
    """在 adata.X 上计算 n_genes / total_counts / pct_counts_mt / pct_counts_ribo。

    为什么在读取阶段计算：确保无论上游是否已做预处理，这些基本 QC 列
    在所有数据集中都已存在且列对齐，避免 stage 2 QC 遇到缺失列分支。

    MT 基因：以 MT- 开头的线粒体基因。
    Ribo 基因：以 RPS/RPL 开头的核糖体蛋白基因。
    """
    # 确保 CSR 格式（scanpy QC 函数的前提）
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    adata.obs["n_genes"] = (
        (adata.X > 0).sum(axis=1).A1
        if sp.issparse(adata.X)
        else (adata.X > 0).sum(axis=1)
    )
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()

    # 线粒体基因比例（MT- 前缀）
    mt_mask = adata.var.index.str.startswith("MT-")
    if mt_mask.any():
        adata.obs["pct_counts_mt"] = (
            np.asarray(adata.X[:, mt_mask].sum(axis=1)).flatten()
            / adata.obs["total_counts"].values
            * 100
        )

    # 核糖体基因比例（RPS/RPL 前缀）
    ribo_mask = adata.var.index.str.startswith(("RPS", "RPL"))
    if ribo_mask.any():
        adata.obs["pct_counts_ribo"] = (
            np.asarray(adata.X[:, ribo_mask].sum(axis=1)).flatten()
            / adata.obs["total_counts"].values
            * 100
        )


# =============================================================================
# manifest 加载与校验
# =============================================================================


def _load_manifest(manifest_path: str) -> dict:
    """从 YAML 文件加载 manifest 并基础校验。"""
    path = Path(manifest_path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest 文件未找到: {manifest_path}")
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    if not isinstance(manifest, dict):
        raise ValueError(f"Manifest {manifest_path} 不是有效的 YAML 映射")
    return manifest


def _validate_manifest(manifest: dict) -> None:
    """校验 manifest schema。致命问题抛出 ValueError。

    校验内容：
    - species（必需，目前仅接受 human）
    - input 块（format + path 必需）
    - source_dataset / project_id / disease_system（必需）
    - original_annotations 段（必需，即使为空列表）
    - qc_overrides 中 skip:true 必须带 reason
    """
    # species 必需，目前仅 human
    species = manifest.get("species")
    if not species:
        raise ValueError("Manifest 缺少必需字段 'species'")
    if species != "human":
        raise ValueError(
            f"物种 '{species}' 不受支持。仅接受 'human'。"
        )

    # input 块必需
    if "input" not in manifest:
        raise ValueError("Manifest 缺少必需段 'input'")

    inp = manifest["input"]
    for key in ("format", "path"):
        if key not in inp:
            raise ValueError(f"Manifest 的 'input' 段缺少必需键 '{key}'")

    fmt = inp["format"]
    supported = {"10x_mtx", "h5ad", "h5", "txt.gz", "rds"}
    if fmt not in supported:
        raise ValueError(
            f"不支持的 input.format '{fmt}'。支持: {', '.join(sorted(supported))}"
        )

    # 数据集标识
    if "source_dataset" not in manifest:
        raise ValueError("Manifest 缺少必需字段 'source_dataset'")
    if "project_id" not in manifest:
        raise ValueError("Manifest 缺少必需字段 'project_id'")
    if "disease_system" not in manifest:
        raise ValueError("Manifest 缺少必需字段 'disease_system'")

    # original_annotations 段必需（即使为空列表）
    if "original_annotations" not in manifest:
        raise ValueError(
            "Manifest 缺少必需段 'original_annotations' "
            "（若无作者标注请写 []）"
        )

    # qc_overrides: skip:true 时 reason 必需
    for step, cfg in manifest.get("qc_overrides", {}).items():
        if cfg.get("skip") and not cfg.get("reason"):
            raise ValueError(
                f"qc_overrides.{step}.skip 为 true 但缺少 'reason'。"
                f"跳过 QC 步骤时必须给出理由。"
            )
