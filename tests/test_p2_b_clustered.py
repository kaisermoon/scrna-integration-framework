"""P2-b-clustered（UX-2）叶子测试：05_clustered.ipynb 的逐 resolution 来源构成、
相邻分辨率稳定性、综合风险标记与研究者显式决策 cell。

测试策略（遵循 notebook 测试纪律）：直接 exec cell 源码，用轻量 mock adata 验证
契约逻辑与决策门禁，不跑真实 scanpy/scVI 重计算。h5ad round-trip 用真实 anndata
但仅针对新增 uns 的列式结构（护住 redline 7：fresh-kernel Run-All 不得因写盘崩溃）。
"""
import ast
import json
import re
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).parents[1]
NOTEBOOK = ROOT / "notebooks/05_clustered.ipynb"


def _source(cell: dict) -> str:
    value = cell.get("source", "")
    return "".join(value) if isinstance(value, list) else value


def _cell(marker: str) -> str:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return next(
        _source(cell) for cell in notebook["cells"]
        if cell["cell_type"] == "code" and marker in _source(cell)
    )


def _base_env(resolutions=None) -> dict:
    """exec PARAMS 拿到全部默认参数，再覆盖 RESOLUTIONS，注入运行时依赖。"""
    resolutions = resolutions if resolutions is not None else [0.2, 0.4]
    env = {
        "pd": pd, "np": np, "sp": sp, "re": re,
        "adjusted_rand_score": adjusted_rand_score,
        "UPSTREAM_RUN_ROOT": "x", "UPSTREAM_RUN_ID": "x",
    }
    exec(_cell("# === PARAMS ==="), env)
    env["RESOLUTIONS"] = resolutions
    return env


def _adata(obs: pd.DataFrame, uns: dict | None = None):
    return SimpleNamespace(obs=obs, uns=dict(uns or {}),
                           n_obs=len(obs), n_vars=2)


# ---------------------------------------------------------------------------
# 静态结构：AST + marker 保留 + 状态标记纪律（P0-i）
# ---------------------------------------------------------------------------
def test_all_code_cells_parse_and_markers_preserved() -> None:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = [_source(c) for c in notebook["cells"] if c["cell_type"] == "code"]
    joined = "\n".join(code)
    ast.parse(joined, filename=str(NOTEBOOK))
    # 改-2 保留 pd.crosstab；shared test 也断言，此处冗余守护。
    assert "pd.crosstab" in joined
    # UX-2 新参数与 uns 键必须落在代码里。
    for token in (
        "BATCH_DOMINANCE_THRESHOLD", "COMPOSITION_KEYS",
        "SELECTED_CLUSTER_KEY", "SELECTION_RATIONALE",
        "05_source_composition", "05_stability_adjacent", "05_resolution_risk",
    ):
        assert token in joined, token


def test_new_cells_use_explicit_completed_markers() -> None:
    # P0-i：新增 cell 用显式 _xxx_completed，不用 dir() 判存在性。
    adj = _cell("# 相邻分辨率稳定性（廉价指标")
    assert adj.splitlines()[0] == "_adjacent_completed = False"
    assert "_adjacent_completed = True" in adj
    assert "in dir()" not in adj


def test_params_defines_ux2_and_decision_defaults() -> None:
    env = _base_env()
    assert env["BATCH_DOMINANCE_THRESHOLD"] == 0.8
    assert env["COMPOSITION_KEYS"] == ["source_dataset", "sample_id", "donor_id", "disease"]
    assert env["SELECTED_CLUSTER_KEY"] is None
    assert env["SELECTION_RATIONALE"] == ""


# ---------------------------------------------------------------------------
# 改-2：逐 resolution 来源构成
# ---------------------------------------------------------------------------
def test_source_composition_flags_single_source_cluster() -> None:
    env = _base_env([0.2])
    # cluster 0 全部来自 ds_a（应标记）；cluster 1 来源均衡（不标记）。
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "0", "0", "1", "1", "1", "1"],
        "source_dataset": ["a", "a", "a", "a", "a", "b", "a", "b"],
    })
    env["adata"] = _adata(obs)
    exec(_cell("# === 跨数据集 cluster 归属偏差检测 ==="), env)
    comp = env["adata"].uns["05_source_composition"]
    assert comp["threshold"] == 0.8
    assert comp["keys_used"] == ["source_dataset"]
    bd = comp["batch_driven_clusters"]
    # 列式结构：等长列表。
    assert set(bd) == {"resolution", "cluster", "dominant_source",
                       "dominant_pct", "n_cells", "flag"}
    assert len({len(v) for v in bd.values()}) == 1
    assert bd["cluster"] == ["0"]
    assert bd["dominant_source"] == ["a"]
    assert bd["flag"] == ["batch_driven"]
    assert bd["n_cells"] == [4]


def test_source_composition_no_source_columns_writes_empty_no_raise() -> None:
    env = _base_env([0.2])
    obs = pd.DataFrame({"leiden_res_0.2": ["0", "0", "1", "1"]})
    env["adata"] = _adata(obs)
    exec(_cell("# === 跨数据集 cluster 归属偏差检测 ==="), env)
    comp = env["adata"].uns["05_source_composition"]
    assert comp["keys_used"] == []
    assert all(len(v) == 0 for v in comp["batch_driven_clusters"].values())


def test_source_composition_prefers_uns_batch_key() -> None:
    env = _base_env([0.2])
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "0", "1"],
        "source_dataset": ["a", "a", "a", "b"],
        "sample_id": ["s1", "s2", "s3", "s4"],
    })
    # 04 写入的 batch_key 指向 source_dataset，应优先作为主判定键。
    env["adata"] = _adata(obs, uns={"harmony_v1": {"batch_key": "source_dataset"}})
    exec(_cell("# === 跨数据集 cluster 归属偏差检测 ==="), env)
    bd = env["adata"].uns["05_source_composition"]["batch_driven_clusters"]
    # cluster 0 全来自 a（>0.8）→ 标记；判定键为 source_dataset。
    assert "0" in bd["cluster"]
    assert set(bd["dominant_source"]) <= {"a", "b"}


# ---------------------------------------------------------------------------
# 新增-A：相邻分辨率稳定性
# ---------------------------------------------------------------------------
def test_adjacent_stability_computes_columnar_uns() -> None:
    env = _base_env([0.2, 0.4, 0.6])
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "1", "1"],
        "leiden_res_0.4": ["0", "1", "2", "2"],
        "leiden_res_0.6": ["0", "1", "2", "3"],
    })
    env["adata"] = _adata(obs)
    exec(_cell("# 相邻分辨率稳定性（廉价指标"), env)
    assert env["_adjacent_completed"] is True
    adj = env["adata"].uns["05_stability_adjacent"]
    assert set(adj) == {"res_low", "res_high", "n_clusters_low",
                        "n_clusters_high", "delta_clusters", "ari"}
    # 3 个分辨率 → 2 个相邻对。
    assert all(len(v) == 2 for v in adj.values())
    assert adj["res_low"] == [0.2, 0.4]
    assert adj["n_clusters_low"] == [2, 3]
    assert adj["n_clusters_high"] == [3, 4]
    assert adj["delta_clusters"] == [1, 1]
    # ARI 为浮点，范围合理。
    assert all(-1.0 <= a <= 1.0 for a in adj["ari"])


# ---------------------------------------------------------------------------
# 新增-B：综合风险标记
# ---------------------------------------------------------------------------
def test_resolution_risk_flags_overfragmentation_and_sample_specific() -> None:
    env = _base_env([0.2, 0.4])
    env["EXPECTED_CLUSTER_MAX"] = 2
    env["MIN_CLUSTER_SIZE"] = 3
    # res 0.4 有 3 个 cluster（> max=2）且含小簇 → over_fragmented。
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "0", "1"],
        "leiden_res_0.4": ["0", "0", "1", "2"],
    })
    # 预置来源构成：res 0.4 有一个 batch_driven cluster。
    uns = {"05_source_composition": {"batch_driven_clusters": {
        "resolution": [0.4], "cluster": ["2"], "dominant_source": ["a"],
        "dominant_pct": [1.0], "n_cells": [1], "flag": ["batch_driven"]}}}
    env["adata"] = _adata(obs, uns=uns)
    exec(_cell("# === 各分辨率综合风险标记"), env)
    risk = env["adata"].uns["05_resolution_risk"]
    assert risk["resolution"] == [0.2, 0.4]
    # res 0.4 过度碎片化，res 0.2 不。
    idx04 = risk["resolution"].index(0.4)
    assert risk["over_fragmented"][idx04] is True
    assert risk["sample_specific_count"][idx04] == 1
    # MARKER_PREVIEW_ENABLED 默认 True → marker 证据可生成。
    assert risk["marker_evidence"][idx04] == "可生成"


def test_resolution_risk_marker_evidence_absent_when_preview_off() -> None:
    env = _base_env([0.2])
    env["MARKER_PREVIEW_ENABLED"] = False
    obs = pd.DataFrame({"leiden_res_0.2": ["0", "0", "1", "1"]})
    env["adata"] = _adata(obs)
    exec(_cell("# === 各分辨率综合风险标记"), env)
    risk = env["adata"].uns["05_resolution_risk"]
    assert risk["marker_evidence"] == ["marker 证据未生成"]


def test_resolution_risk_tolerates_missing_source_composition() -> None:
    # 顺序解耦：来源构成 cell 未运行时用 .get 容错，不报错。
    env = _base_env([0.2])
    obs = pd.DataFrame({"leiden_res_0.2": ["0", "0", "1", "1"]})
    env["adata"] = _adata(obs)  # uns 无 05_source_composition
    exec(_cell("# === 各分辨率综合风险标记"), env)
    risk = env["adata"].uns["05_resolution_risk"]
    assert risk["sample_specific_count"] == [0]


# ---------------------------------------------------------------------------
# 新增-C：研究者显式决策 cell
# ---------------------------------------------------------------------------
def _decision_env(selected, resolutions=None, obs=None):
    env = _base_env(resolutions or [0.2, 0.4])
    env["SELECTED_CLUSTER_KEY"] = selected
    env["SELECTION_RATIONALE"] = "res 0.4 匹配预期大类数"
    if obs is None:
        obs = pd.DataFrame({
            "leiden_res_0.2": ["0", "0", "1", "1"],
            "leiden_res_0.4": ["0", "1", "2", "2"],
        })
    env["adata"] = _adata(obs)
    return env


def test_decision_none_keeps_needs_review_no_raise() -> None:
    env = _decision_env(None)
    exec(_cell("# === 研究者显式选择 selected_cluster_key"), env)
    # 未选择：不写 selected_cluster_key，不 raise（护 Run-All / CI）。
    assert "selected_cluster_key" not in env["adata"].uns


def test_decision_valid_writes_uns() -> None:
    env = _decision_env("leiden_res_0.4")
    exec(_cell("# === 研究者显式选择 selected_cluster_key"), env)
    assert env["adata"].uns["selected_cluster_key"] == "leiden_res_0.4"
    assert env["adata"].uns["selection_rationale"] == "res 0.4 匹配预期大类数"


def test_decision_key_not_in_obs_raises() -> None:
    env = _decision_env("leiden_res_9.9")
    with pytest.raises(ValueError, match="not in obs"):
        exec(_cell("# === 研究者显式选择 selected_cluster_key"), env)


def test_decision_resolution_not_swept_raises() -> None:
    # 列存在于 obs，但其 resolution 不在 RESOLUTIONS 中。
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "1", "1"],
        "leiden_res_0.4": ["0", "1", "2", "2"],
        "leiden_res_1.5": ["0", "1", "2", "3"],
    })
    env = _decision_env("leiden_res_1.5", resolutions=[0.2, 0.4], obs=obs)
    with pytest.raises(ValueError, match="not swept"):
        exec(_cell("# === 研究者显式选择 selected_cluster_key"), env)


def test_decision_single_cluster_raises() -> None:
    obs = pd.DataFrame({
        "leiden_res_0.2": ["0", "0", "0", "0"],
        "leiden_res_0.4": ["0", "1", "2", "2"],
    })
    env = _decision_env("leiden_res_0.2", obs=obs)
    with pytest.raises(ValueError):
        exec(_cell("# === 研究者显式选择 selected_cluster_key"), env)


def test_decision_does_not_auto_assign_selected_key() -> None:
    # 决策 cell 绝不根据指标自动给 SELECTED_CLUSTER_KEY 赋计算值。
    source = _cell("# === 研究者显式选择 selected_cluster_key")
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                assert not (isinstance(tgt, ast.Name)
                            and tgt.id == "SELECTED_CLUSTER_KEY")


# ---------------------------------------------------------------------------
# 改-3：manifest 追加键
# ---------------------------------------------------------------------------
def test_manifest_appends_selection_and_composition_keys() -> None:
    final = _cell("# === Stage 05 draft checkpoint")
    assert '"selected_cluster_key": SELECTED_CLUSTER_KEY' in final
    assert '"selection_rationale": SELECTION_RATIONALE' in final
    assert '"source_composition_summary"' in final
    # 追加不得引入 promote_run（P0/P1 红线）。
    assert "promote_run" not in final


# ---------------------------------------------------------------------------
# redline 7 守护：新增 uns 列式结构必须能安全 h5ad round-trip
# ---------------------------------------------------------------------------
def test_new_uns_columnar_dicts_survive_h5ad_roundtrip(tmp_path: Path) -> None:
    adata = ad.AnnData(sp.csr_matrix(np.ones((4, 2), dtype=np.float32)))
    adata.uns["05_source_composition"] = {
        "threshold": 0.8, "keys_used": ["source_dataset"],
        "batch_driven_clusters": {
            "resolution": [0.4], "cluster": ["2"], "dominant_source": ["a"],
            "dominant_pct": [0.9], "n_cells": [10], "flag": ["batch_driven"]},
    }
    adata.uns["05_stability_adjacent"] = {
        "res_low": [0.2], "res_high": [0.4], "n_clusters_low": [2],
        "n_clusters_high": [3], "delta_clusters": [1], "ari": [0.8]}
    adata.uns["05_resolution_risk"] = {
        "resolution": [0.2, 0.4], "n_clusters": [2, 3], "n_small_clusters": [0, 1],
        "over_fragmented": [False, True], "sample_specific_count": [0, 1],
        "marker_evidence": ["可生成", "可生成"], "doublet_enriched": [False, False]}
    path = tmp_path / "roundtrip.h5ad"
    adata.write_h5ad(path)  # 不得抛 TypeError（list-of-dict 会抛）
    restored = ad.read_h5ad(path)
    assert "05_source_composition" in restored.uns
    assert "05_stability_adjacent" in restored.uns
    assert "05_resolution_risk" in restored.uns
