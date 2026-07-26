"""IO 模块：基因 ID 双向同步、基因组位置注入、batch 诊断。

src/notebook 边界铁律（2026-07-10 PI 定稿）：进 src/ 的仅为技术管道。
本模块保留三组技术管道函数：
- sync_gene_ids + 三个私有 helper：基因 ID symbol↔ensembl 双向同步（mygene 在线查询）
- inject_genomic_positions：向 adata.var 注入 chromosome/start/end（CNV 推断等下游分析用）
- summarize_batch_keys：batch 键诊断打印（确定性，零 LLM）

面向非计算机专业 PI/学生：从上到下线性可读，打开文件就能看懂每个函数的用途。
"""

from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd

# =============================================================================
# 基因 ID 双向同步（mygene 在线查询 + 批量处理）
# =============================================================================


def _is_ensembl_id(gene_id: str) -> bool:
    """判断 gene_id 是否为 Ensembl ID 格式。

    Ensembl ID 以 "ENS" 开头（ENSG/ENSMUSG/ENSRNOG/ENSDARG 等），
    兼容人类、小鼠、大鼠、斑马鱼等常见模式生物。
    大小写不敏感（转大写后判断）。
    """
    return str(gene_id).upper().startswith("ENS")


def _infer_species_from_ensembl_prefix(gene_id: str) -> str:
    """从单个 Ensembl ID 前缀推断 mygene species 参数。

    Returns
    -------
    str
        "human" / "mouse" / "rat" / "all"（无法判断时返回 "all"，
        mygene 在此模式下不加 species 过滤，通配所有物种）。
    """
    gid = str(gene_id).upper()
    if gid.startswith("ENSMUSG"):
        return "mouse"
    if gid.startswith("ENSRNOG"):
        return "rat"
    if gid.startswith("ENSG"):
        return "human"
    return "all"


def _infer_species_from_batch_ensembl_ids(gene_ids, adata=None) -> str:
    """从一批 Ensembl ID 推断多数（majority）物种，用于 mygene 批量查询。

    采样前 20 个基因 ID，按前缀统计物种，返回多数派。
    若无法判断（全部为 non-ENS 或出现平票），回退到 adata.uns['species']
    再做归一化映射，最终回退 'human'。

    Parameters
    ----------
    gene_ids : iterable
        基因 ID 列表（可以是 Ensembl ID 或 symbol）。
    adata : AnnData, optional
        当多数投票失败时，从此对象的 adata.uns['species'] 读取物种信息。

    Returns
    -------
    str
        mygene querymany(species=...) 参数值。
    """
    from collections import Counter

    species_votes = []
    for gid in list(gene_ids)[:20]:
        s = _infer_species_from_ensembl_prefix(gid)
        if s != "all":
            species_votes.append(s)

    if species_votes:
        return Counter(species_votes).most_common(1)[0][0]

    # 多数投票失败 → 回退 adata.uns['species']（做归一化映射）
    if adata is not None and "species" in adata.uns:
        sp_raw = str(adata.uns["species"]).lower()
        _species_norm = {
            "human": "human", "homo_sapiens": "human",
            "mouse": "mouse", "mus_musculus": "mouse",
            "rat": "rat", "rattus_norvegicus": "rat",
        }
        return _species_norm.get(sp_raw, "human")

    # 最终回退 human
    return "human"


def sync_gene_ids(adata, gene_id_format: str = "auto") -> None:
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
    is_ensembl = _is_ensembl_id(idx_sample)

    if gene_id_format == "symbol" or (gene_id_format == "auto" and not is_ensembl):
        _sync_symbol_to_ensembl(adata)
    elif gene_id_format == "ensembl" or (gene_id_format == "auto" and is_ensembl):
        _sync_ensembl_to_symbol(adata)
    else:
        _sync_symbol_to_ensembl(adata)


def _sync_symbol_to_ensembl(adata) -> None:
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
            gene_ids.astype(str).str.upper().str.startswith("ENS"), ""
        )
        n_mapped = (adata.var["ensembl_id"] != "").sum()
        if n_mapped > 0:
            return

    # 回退：在线查询 mygene
    _mygene_symbol_to_ensembl(adata)


def _mygene_symbol_to_ensembl(adata) -> None:
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


def _sync_ensembl_to_symbol(adata) -> None:
    """var.index 是 Ensembl ID；转换为 gene symbol。

    三层优先级策略：
    ① feature_name 列（CELLxGENE 惯例，离线可用，主路径）
    ② 对 feature_name 中仍未转换的 Ensembl ID（lncRNA/假基因），
       用 mygene 在线查询二次补全
    ③ 无 feature_name 列时 mygene 全量兜底在线查询

    每次设置 var_names 后紧跟去重（scanpy 生态硬要求 var.index 唯一），
    确保同一 gene symbol 对应多条不同 Ensembl ID 时不会因重复索引报错。
    """
    # 幂等检查：ensembl_id 列已存在且非全空 → 跳过重复转换
    if "ensembl_id" in adata.var.columns and adata.var["ensembl_id"].notna().any():
        return

    # —— ① feature_name 列（CELLxGENE 惯例，离线，主路径） ——
    # CELLxGENE 导出的 h5ad 通常有 feature_name 列存储 gene symbol，
    # 优先使用它——离线即可用，无需网络，且已处理无 symbol 基因的 Ensembl ID 回退。
    if "feature_name" in adata.var.columns:
        feature_names = adata.var["feature_name"].values.astype(str)
        n_total = len(feature_names)
        n_has_symbol = sum(1 for v in feature_names if not _is_ensembl_id(v))

        if n_has_symbol > 0:
            # 保存原始 Ensembl ID 到 ensembl_id 列
            original_ensembl = adata.var.index.values.astype(str)

            # 用 feature_name 的值替换 var.index（无 symbol 的基因保留 Ensembl ID）
            adata.var_names = feature_names

            # 处理重复 gene symbol：对重名追加 -1、-2 等后缀，
            # 确保 var.index 唯一（scanpy 生态硬要求）
            if not adata.var.index.is_unique:
                n_dup = int(adata.var.index.duplicated().sum())
                adata.var_names_make_unique()
                print(
                    f"基因 ID 同步：{n_dup} 个重复 gene symbol 已自动去重 "
                    f"（追加 -N 后缀）"
                )

            adata.var["ensembl_id"] = original_ensembl

            n_retained = n_total - n_has_symbol
            print(
                f"基因 ID 同步（feature_name）: "
                f"{n_has_symbol}/{n_total} 个基因转换为 gene symbol，"
                f"{n_retained} 个基因无已知 symbol 保留 Ensembl ID"
            )

            # —— ② mygene 二次补全：对 feature_name 中仍是 Ensembl ID 的基因 ——
            # 为什么需要这一步：CELLxGENE 的 feature_name 对 lncRNA/假基因等
            # 可能仍保留 Ensembl ID，mygene 在线查询可补全这部分基因的 symbol。
            _ensg_indices = [
                i for i, g in enumerate(adata.var.index)
                if _is_ensembl_id(str(g))
            ]
            _n_ensg_remaining = len(_ensg_indices)

            if _n_ensg_remaining > 0:
                import mygene

                mg = mygene.MyGeneInfo()
                # 去掉 var_names_make_unique() 追加的 -N 后缀后查询，
                # 否则拼接后的 ID（如 "ENSG00000123456-1"）在 mygene 中查不到
                _ensg_to_query = [
                    re.sub(r"-\d+$", "", str(adata.var.index[i]))
                    for i in _ensg_indices
                ]
                symbol_map: dict[str, str] = {}
                _batch_species = _infer_species_from_batch_ensembl_ids(
                    _ensg_to_query, adata
                )

                batch_size = 1000
                for batch_start in range(0, len(_ensg_to_query), batch_size):
                    batch = _ensg_to_query[
                        batch_start: batch_start + batch_size
                    ]
                    try:
                        results = mg.querymany(
                            batch, scopes="ensembl.gene", fields="symbol",
                            species=_batch_species, returnall=False,
                        )
                    except Exception:
                        warnings.warn(
                            f"mygene 补全查询失败（批次 {batch_start}）；"
                            f"保留 Ensembl ID 作为索引",
                            stacklevel=2,
                        )
                        continue
                    for r in results:
                        query = r.get("query", "")
                        sym = r.get("symbol")
                        if sym:
                            symbol_map[query] = sym

                if symbol_map:
                    new_index = list(adata.var.index)
                    n_mapped = 0
                    for i in _ensg_indices:
                        clean = re.sub(
                            r"-\d+$", "", str(adata.var.index[i])
                        )
                        if clean in symbol_map:
                            new_index[i] = symbol_map[clean]
                            n_mapped += 1

                    adata.var_names = new_index

                    # mygene 补全后再次检查重复：
                    # 新返回的 symbol 可能与已有 symbol 冲突
                    if not adata.var.index.is_unique:
                        n_dup_after = int(
                            adata.var.index.duplicated().sum()
                        )
                        adata.var_names_make_unique()
                        print(
                            f"mygene 补全后去重：{n_dup_after} 个重复 symbol "
                            f"追加 -N 后缀"
                        )

                    print(
                        f"mygene 补全: {n_mapped} 个基因新增 symbol"
                        f"（共 {_n_ensg_remaining} 个 ENSG 基因查询）"
                    )
                else:
                    print(
                        f"mygene 补全: 0 个基因新增 symbol"
                        f"（共 {_n_ensg_remaining} 个 ENSG 基因查询，"
                        f"全部未命中）"
                    )

            return
        # feature_name 全是 Ensembl 格式（无真实 symbol）→ 回退 mygene 全量兜底
        # 注意：此时尚未修改 var.index 或添加 ensembl_id 列，无需回滚

    # —— ③ 无 feature_name 列时 mygene 全量兜底 ——
    import mygene

    mg = mygene.MyGeneInfo()
    ensembl_ids = list(adata.var.index)
    symbol_map_result: dict[str, str] = {}
    _fallback_species = _infer_species_from_batch_ensembl_ids(
        ensembl_ids, adata
    )

    batch_size = 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i: i + batch_size]
        try:
            results = mg.querymany(
                batch, scopes="ensembl.gene", fields="symbol",
                species=_fallback_species, returnall=False,
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
    adata.var_names = new_index

    # mygene 查询后检查重复：同一 symbol 可能对应多个 Ensembl ID
    if not adata.var.index.is_unique:
        n_dup = int(adata.var.index.duplicated().sum())
        adata.var_names_make_unique()
        print(
            f"基因 ID 同步（mygene）: {n_dup} 个重复 gene symbol "
            f"已自动去重（追加 -N 后缀）"
        )

    adata.var["ensembl_id"] = adata.var["_orig_ensembl"]
    del adata.var["_orig_ensembl"]

    n_missing = sum(1 for v in new_index if _is_ensembl_id(str(v)))
    if n_missing > 0:
        warnings.warn(
            f"基因 ID 同步: {n_missing}/{len(ensembl_ids)} 个 Ensembl ID 无法通过 "
            f"mygene 映射到 gene symbol。保留 Ensembl ID 作为索引，"
            f"对应基因的 ensembl_id 留空。",
            stacklevel=2,
        )
        mask = [_is_ensembl_id(str(v)) for v in adata.var.index]
        adata.var.loc[mask, "ensembl_id"] = ""


# =============================================================================
# 基因组位置注入（公开，用于 CNV 推断等下游分析）
# =============================================================================


def inject_genomic_positions(
    adata,
    species: str = "human",
    batch_size: int = 1000,
):
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
# batch 键诊断（公开，纯确定性打印）
# =============================================================================


def summarize_batch_keys(adata) -> None:
    """打印合并后 AnnData 各 batch 相关列的取值分布摘要。

    供 PI 在多源合并后手动调用，排查 batch 键语义不一致等问题。
    纯确定性——零 LLM。不改动数据，仅打印诊断信息。

    为什么需要独立函数：batch 语义校验
    发生在 anndata.concat 之后（per-dataset notebook 合并 cell 里），不在单源
    读取内。合并后的 adata 才具备完整的 batch 列对比上下文。

    Parameters
    ----------
    adata : AnnData
        合并后的 AnnData 对象（如 stage1 多源合并后）。
    """
    batch_cols = ["batch", "sample_id", "source_dataset"]
    for col in batch_cols:
        if col not in adata.obs.columns:
            print(f"[summarize_batch_keys] 列 '{col}' 不存在")
            continue
        series = adata.obs[col]
        n_unique = series.nunique()
        has_nan = bool(series.isna().any())
        dtype_str = str(series.dtype)
        print(
            f"[summarize_batch_keys] {col}: dtype={dtype_str}, "
            f"n_unique={n_unique}, has_NaN={has_nan}"
        )
        # 整数批次号与长字符串 ID 混存提示
        if n_unique > 1:
            sample_vals = series.dropna().unique()[:10]
            str_vals = [v for v in sample_vals if not str(v).isdigit()]
            int_vals = [v for v in sample_vals if str(v).isdigit()]
            if str_vals and int_vals:
                print(
                    f"  ⚠ 整数与字符串混存："
                    f"整数示例={[str(v) for v in int_vals[:3]]}, "
                    f"字符串示例={[str(v) for v in str_vals[:3]]}"
                )
