"""UX-3 run 管理地基：run_contract.py 新增 helper 的叶子测试。

测试策略：
- 所有 manifest 夹具先过 validate_manifest 再落盘（atomic_write_json），
  确保枚举器读的是真实 v1 manifest 而非 stub。
- 不用 git show <sha>/git 依赖，tmp_path 自足，CI shallow-clone 安全。
- 每个测试族含非空跑证明断言（若被测函数返回常量/占位，则断言必红）。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from scrna_integration.run_contract import (
    LARGE_MATRIX_SUFFIXES,
    MANIFEST_SCHEMA_VERSION,
    PINNED_MARKER_NAME,
    RunCategory,
    RunManifestRecord,
    StageStatus,
    atomic_write_json,
    classify_run,
    diff_effective_parameters,
    enumerate_cleanup_candidates,
    enumerate_run_manifests,
    sha256_file,
    validate_manifest,
)

# ---------------------------------------------------------------------------
# v1 manifest 工厂：产出通过 validate_manifest 的完整 manifest dict
# ---------------------------------------------------------------------------

def _v1_manifest(
    *,
    run_id: str = "run-001",
    stage: str = "04_embedded",
    status: str = "SUCCESS",
    effective_parameters: dict | None = None,
    checkpoint_path: str = "checkpoint.h5ad",
    checkpoint_sha256: str = "0" * 64,
    failed: bool = False,
) -> dict:
    """构建通过 validate_manifest(mapping, state="draft") 的完整 v1 manifest dict。

    调用方落盘后须确保 checkpoint_path 对应的文件真实存在且 hash 匹配，
    否则 validate_manifest 文件模式会报 hash 不匹配。
    """
    manifest: dict = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "run_id": run_id,
        "stage": stage,
        "stage_status": status,
        "started_at": "2026-07-14T01:02:03Z",
        "completed_at": "2026-07-14T01:02:04.123Z",
        "inputs": [
            {"role": "upstream", "path": "upstream/manifest.json",
             "kind": "file", "sha256": "1" * 64},
        ],
        "effective_parameters": effective_parameters if effective_parameters is not None else {},
        "runtime_provenance": {
            "python": "3.11",
            "platform": "test",
            "git_available": False,
            "git_commit": "unavailable",
            "git_dirty": None,
            "git_status_sha256": None,
            "git_untracked_count": None,
            "git_tracked_dirty": None,
            "git_tracked_diff_sha256": None,
            "packages": {},
        },
        "method_status": {"pca": "success"},
        "hard_postconditions": {"checkpoint_written": True},
        "warnings": ["warning"] if status == "SUCCESS_WITH_WARNINGS" else [],
        "artifacts": [],
    }
    if failed:
        manifest["failure"] = {"type": "RuntimeError", "message": "simulated failure"}
    else:
        manifest["checkpoint"] = {"path": checkpoint_path, "sha256": checkpoint_sha256}
    return manifest


def _write_valid_draft(
    run_dir: Path,
    manifest_dict: dict,
    *,
    checkpoint_content: bytes = b"checkpoint",
    state: str = "draft",
) -> Path:
    """将 manifest dict 写入 state_dir（draft 或 promoted），并写入真实 checkpoint 文件。

    先用 validate_manifest(mapping, state=...) 校验 manifest_dict，
    确保夹具是真实 v1 manifest 而非 stub，再 atomic_write_json 落盘。
    checkpoint 文件按 manifest["checkpoint"]["path"] 创建，hash 自动匹配。

    promoted 场景：先写入 draft 目录，再通过 os.replace 模拟原子提升，
    使其符合 promote_run 的目录布局（draft → promoted 整目录重命名）。
    """
    # 深拷贝：防止 validate_manifest 修改传入的 dict
    manifest = json.loads(json.dumps(manifest_dict))

    # 先通过 mapping 模式 schema 校验
    validate_manifest(manifest_dict, state=state)

    # 确保 run_id 匹配目录名
    manifest["run_id"] = run_dir.name

    draft_dir = run_dir / "draft"
    promoted_dir = run_dir / "promoted"

    # 写入 draft（promoted 也先写 draft，最后再 rename）
    draft_dir.mkdir(parents=True, exist_ok=True)

    # 写 checkpoint 文件
    if "checkpoint" in manifest:
        cp_path = draft_dir / manifest["checkpoint"]["path"]
        cp_path.parent.mkdir(parents=True, exist_ok=True)
        cp_path.write_bytes(checkpoint_content)
        manifest["checkpoint"]["sha256"] = sha256_file(cp_path)

    # 写入 manifest 到 draft
    atomic_write_json(draft_dir / "manifest.json", manifest)

    if state == "promoted":
        # simulate promote_run: os.replace(draft, promoted)
        os.replace(draft_dir, promoted_dir)
        return promoted_dir
    else:
        return draft_dir


# ---------------------------------------------------------------------------
# 便捷 fixture 别名（从 RunCategory 映射，减少魔法字符串）
# ---------------------------------------------------------------------------

# 直接引用 RunCategory 值作为模块常量
FAILED = RunCategory.FAILED
PROMOTED = RunCategory.PROMOTED
PINNED = RunCategory.PINNED
SUPERSEDED = RunCategory.SUPERSEDED
SUCCESS = StageStatus.SUCCESS
SUCCESS_WITH_WARNINGS = StageStatus.SUCCESS_WITH_WARNINGS
NEEDS_REVIEW = StageStatus.NEEDS_REVIEW


# ===================================================================
# T1  classify_run — 纯映射
# ===================================================================

class TestClassifyRun:
    """classify_run 的四态映射：FAILED > promoted > pinned > superseded。"""

    @pytest.mark.parametrize(
        ("stage_status", "promoted", "pinned", "expected"),
        [
            ("FAILED", False, False, FAILED),
            ("SUCCESS", True, False, PROMOTED),
            ("SUCCESS", False, True, PINNED),
            ("SUCCESS", False, False, SUPERSEDED),
            # 活跃 draft（NEEDS_REVIEW）落入 superseded 桶，细粒度由 stage_status 承载
            ("NEEDS_REVIEW", False, False, SUPERSEDED),
            # 优先级验证：FAILED 覆盖 pinned
            ("FAILED", False, True, FAILED),
            # 优先级验证：promoted 覆盖 pinned
            ("SUCCESS", True, True, PROMOTED),
        ],
    )
    def test_classification(self, stage_status, promoted, pinned, expected):
        assert classify_run(stage_status, promoted=promoted, pinned=pinned) == expected

    def test_invalid_stage_status_raises_valueerror(self):
        with pytest.raises(ValueError, match="invalid stage_status"):
            classify_run("BOGUS", promoted=False, pinned=False)

    def test_non_empty_proof_distinct_branches_return_different_values(self):
        """非空跑证明：四个分支必须返回互不相同的值。

        若 classify_run 被写成返回常量，此断言必红。
        """
        r1 = classify_run("FAILED", promoted=False, pinned=False)
        r2 = classify_run("SUCCESS", promoted=True, pinned=False)
        r3 = classify_run("SUCCESS", promoted=False, pinned=True)
        r4 = classify_run("SUCCESS", promoted=False, pinned=False)
        results = {r1, r2, r3, r4}
        assert len(results) == 4, (
            f"四个不同分支应返回四种不同 RunCategory，实际得到 {len(results)} 种: {results}"
        )


# ===================================================================
# T2  enumerate_run_manifests — 多态 RUN_ROOT
# ===================================================================

class TestEnumerateRunManifests:
    """enumerate_run_manifests 枚举多态 run 目录。

    非空跑证明：构造四个具有不同 stage_status/category 的 run，
    逐一断言各 record 的具体字段值——若 enumerate 返回常量/占位，必红。
    """

    def _build_root(self, tmp_path: Path) -> Path:
        """构建一个含四态 run + 干扰目录的 RUN_ROOT。

        布局：
          run_A/  promoted + SUCCESS（有真实 checkpoint 文件）
          run_B/  draft + SUCCESS + pin marker
          run_C/  draft + FAILED（failure 字段，无 checkpoint）
          run_D/  draft + NEEDS_REVIEW
          not_a_run/  无 manifest 的子目录（干扰）
          .hidden/  点开头的隐藏目录（干扰）
        """
        root = tmp_path / "runs"
        root.mkdir()

        # run_A: promoted + SUCCESS + checkpoint (3 字节)
        run_a_dir = root / "run_A"
        run_a_dir.mkdir()
        ma = _v1_manifest(
            run_id="run_A", status="SUCCESS",
            effective_parameters={"BATCH_KEY": "sample", "N_HVG": 2000},
        )
        _write_valid_draft(run_a_dir, ma, checkpoint_content=b"ABC", state="promoted")
        # run_B: draft + SUCCESS + pin marker
        run_b_dir = root / "run_B"
        run_b_dir.mkdir()
        mb = _v1_manifest(
            run_id="run_B", status="SUCCESS",
            effective_parameters={"BATCH_KEY": "sample", "N_HVG": 3000},
        )
        _write_valid_draft(run_b_dir, mb, checkpoint_content=b"DEFGH")

        # 建 pin marker
        (run_b_dir / PINNED_MARKER_NAME).write_text("")

        # run_C: draft + FAILED（无 checkpoint）
        run_c_dir = root / "run_C"
        run_c_dir.mkdir()
        mc = _v1_manifest(run_id="run_C", status="FAILED", failed=True)
        mc.pop("checkpoint", None)
        # 先校验 mapping
        validate_manifest(mc, state="draft")
        draft_c = run_c_dir / "draft"
        draft_c.mkdir(parents=True)
        atomic_write_json(draft_c / "manifest.json", mc)

        # run_D: draft + NEEDS_REVIEW
        run_d_dir = root / "run_D"
        run_d_dir.mkdir()
        md = _v1_manifest(
            run_id="run_D", status="NEEDS_REVIEW",
            effective_parameters={"BATCH_KEY": "sample"},
        )
        _write_valid_draft(run_d_dir, md, checkpoint_content=b"XYZ")

        # 干扰目录：不是 run
        (root / "not_a_run").mkdir()
        (root / ".hidden").mkdir()

        return root

    @pytest.fixture
    def run_root(self, tmp_path: Path) -> Path:
        return self._build_root(tmp_path)

    def test_finds_all_four_runs(self, run_root):
        records = enumerate_run_manifests(run_root)
        run_ids = {r.run_id for r in records}
        assert run_ids == {"run_A", "run_B", "run_C", "run_D"}

    def test_run_a_promoted_with_checkpoint(self, run_root):
        records = {r.run_id: r for r in enumerate_run_manifests(run_root)}
        rec = records["run_A"]
        assert rec.state == "promoted"
        assert rec.category == PROMOTED
        assert rec.stage_status == "SUCCESS"
        assert rec.pinned is False
        assert rec.effective_parameters == {"BATCH_KEY": "sample", "N_HVG": 2000}
        assert rec.error is None
        # 输出文件：已知字节数 3（b"ABC"）
        assert rec.output_path is not None
        assert rec.output_size_bytes == 3, (
            f"output_size_bytes 应来自 os.stat 真读（3），实际 {rec.output_size_bytes}"
        )
        assert rec.output_mtime is not None

    def test_run_b_draft_pinned(self, run_root):
        records = {r.run_id: r for r in enumerate_run_manifests(run_root)}
        rec = records["run_B"]
        assert rec.state == "draft"
        assert rec.category == PINNED
        assert rec.pinned is True
        assert rec.stage_status == "SUCCESS"
        assert rec.output_path is not None
        assert rec.output_size_bytes == 5  # b"DEFGH"

    def test_run_c_failed_no_checkpoint(self, run_root):
        records = {r.run_id: r for r in enumerate_run_manifests(run_root)}
        rec = records["run_C"]
        assert rec.state == "draft"
        assert rec.category == FAILED
        assert rec.stage_status == "FAILED"
        # FAILED run：output_path/size/mtime 全部 None
        assert rec.output_path is None
        assert rec.output_size_bytes is None
        assert rec.output_mtime is None

    def test_run_d_needs_review_superseded(self, run_root):
        records = {r.run_id: r for r in enumerate_run_manifests(run_root)}
        rec = records["run_D"]
        assert rec.state == "draft"
        assert rec.category == SUPERSEDED
        # 细粒度状态保留
        assert rec.stage_status == "NEEDS_REVIEW"

    def test_skips_non_run_directories(self, run_root):
        records = enumerate_run_manifests(run_root)
        run_ids = {r.run_id for r in records}
        assert "not_a_run" not in run_ids
        # 隐藏目录也不应是 run
        assert not any(r.run_id.startswith(".") for r in records)

    def test_bad_manifest_produces_error_record(self, tmp_path):
        root = tmp_path / "runs"
        root.mkdir()
        run_dir = root / "bad-run"
        run_dir.mkdir()
        draft_dir = run_dir / "draft"
        draft_dir.mkdir(parents=True)
        # 写入非法 JSON
        (draft_dir / "manifest.json").write_text("this is not json", encoding="utf-8")

        records = enumerate_run_manifests(root)
        assert len(records) == 1
        assert records[0].run_id == "bad-run"
        assert records[0].error is not None
        assert "parse error" in records[0].error

    def test_non_empty_proof_four_categories_distinct(self, run_root):
        """非空跑证明：四个 run 的 category 必须各不相同。

        若 enumerate 返回常量 category，此断言必红。
        """
        records = enumerate_run_manifests(run_root)
        categories = {r.run_id: r.category for r in records}
        assert categories["run_A"] == PROMOTED
        assert categories["run_B"] == PINNED
        assert categories["run_C"] == FAILED
        assert categories["run_D"] == SUPERSEDED
        # 额外保障：四个值互不相同
        cat_set = set(categories.values())
        assert len(cat_set) == 4, f"四个 run 应有四种不同 category，实际: {cat_set}"


# ===================================================================
# T3  diff_effective_parameters — 跨 run 参数结构化 diff
# ===================================================================

class TestDiffEffectiveParameters:
    """diff_effective_parameters 的跨 run 比较，覆盖值差异/缺失/嵌套/边界。"""

    def test_differing_value_across_two_runs(self):
        result = diff_effective_parameters({
            "run_A": {"BATCH_KEY": "sample", "N_HVG": 2000},
            "run_B": {"BATCH_KEY": "sample", "N_HVG": 3000},
        })
        assert result["run_ids"] == ["run_A", "run_B"]
        assert "N_HVG" in result["differing_keys"]
        assert "BATCH_KEY" not in result["differing_keys"]
        assert "BATCH_KEY" in result["shared_keys"]
        # N_HVG 在两 run 中 present=True，值不同
        n_hvg = result["parameters"]["N_HVG"]
        assert n_hvg["differs"] is True
        assert n_hvg["values"]["run_A"] == {"present": True, "value": 2000}
        assert n_hvg["values"]["run_B"] == {"present": True, "value": 3000}

    def test_key_present_in_one_run_only(self):
        result = diff_effective_parameters({
            "run_A": {"EXTRA_KEY": 42},
            "run_B": {},
        })
        assert "EXTRA_KEY" in result["differing_keys"]
        extra = result["parameters"]["EXTRA_KEY"]
        assert extra["differs"] is True
        assert extra["values"]["run_A"] == {"present": True, "value": 42}
        assert extra["values"]["run_B"] == {"present": False}

    def test_nested_dict_deep_compare(self):
        """嵌套 dict 深比较：内层值不同应识别为 differ。"""
        result = diff_effective_parameters({
            "run_A": {"MODEL_PARAMS": {"lr": 0.001, "epochs": 100}},
            "run_B": {"MODEL_PARAMS": {"lr": 0.001, "epochs": 200}},
        })
        assert "MODEL_PARAMS" in result["differing_keys"]
        # 相同嵌套 → 不 differ
        result2 = diff_effective_parameters({
            "run_A": {"MODEL_PARAMS": {"lr": 0.001}},
            "run_B": {"MODEL_PARAMS": {"lr": 0.001}},
        })
        assert "MODEL_PARAMS" in result2["shared_keys"]

    def test_empty_input_returns_empty_structure(self):
        result = diff_effective_parameters({})
        assert result == {"run_ids": [], "parameters": {}, "differing_keys": [], "shared_keys": []}

    def test_single_run_has_no_diff(self):
        result = diff_effective_parameters({"run_A": {"KEY": 1, "OTHER": 2}})
        assert result["differing_keys"] == []
        assert set(result["shared_keys"]) == {"KEY", "OTHER"}

    def test_non_empty_proof_same_params_no_diff(self):
        """非空跑证明：完全相同的参数集 → differing_keys 为空。"""
        params = {"A": 1, "B": {"nested": True}}
        result = diff_effective_parameters({"r1": params, "r2": params})
        assert result["differing_keys"] == []

    def test_non_empty_proof_one_change_detected(self):
        """非空跑证明：仅改一个值 → differing_keys 精确命中所改 key。"""
        result = diff_effective_parameters({
            "r1": {"X": 1, "Y": 2},
            "r2": {"X": 1, "Y": 99},
        })
        assert result["differing_keys"] == ["Y"], (
            f"只改了 Y，differing_keys 应仅含 Y，实际: {result['differing_keys']}"
        )


# ===================================================================
# T4  enumerate_cleanup_candidates — 清理候选枚举（只读）
# ===================================================================

class TestEnumerateCleanupCandidates:
    """enumerate_cleanup_candidates：只枚举 superseded/failed 的大文件，不删除。

    非空跑证明：promoted/pinned 的大 h5ad 不出现在结果；
    调用后所有候选文件仍存在（证明未误删）。
    """

    def _make_record(
        self,
        run_id: str,
        category: RunCategory,
        state_dir: Path,
    ) -> RunManifestRecord:
        return RunManifestRecord(
            run_id=run_id,
            state="draft",
            stage="04_embedded",
            stage_status="SUCCESS",
            category=category,
            pinned=(category == PINNED),
            effective_parameters={},
            manifest_path=state_dir / "manifest.json",
            state_dir=state_dir,
            output_path=None,
            output_size_bytes=None,
            output_mtime=None,
            error=None,
        )

    def _make_h5ad(self, parent: Path, name: str, size: int) -> Path:
        path = parent / name
        parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"x" * size)
        return path

    @pytest.fixture
    def record_set(self, tmp_path: Path) -> list[RunManifestRecord]:
        """构造四个 run，各有 h5ad 文件在不同的 state_dir。"""
        root = tmp_path / "records"
        root.mkdir()

        records = []
        for rid, cat in [
            ("promoted-1", PROMOTED),
            ("pinned-1", PINNED),
            ("superseded-1", SUPERSEDED),
            ("failed-1", FAILED),
        ]:
            sd = root / rid / "draft"
            sd.mkdir(parents=True)
            rec = self._make_record(rid, cat, sd)
            # 每个 record 放一个大 h5ad
            self._make_h5ad(sd, "output.h5ad", 1024)
            records.append(rec)

        # superseded run 中额外放一个嵌套子目录的 h5ad（测试递归扫描）
        extra_dir = root / "superseded-1" / "draft" / "nested"
        extra_dir.mkdir(parents=True)
        self._make_h5ad(extra_dir, "deep.h5ad", 2048)

        # failed run 中放 .txt（不应被收集，因为后缀不匹配）
        self._make_h5ad(root / "failed-1" / "draft", "log.txt", 9999)

        return records

    def test_only_superseded_and_failed_included(self, record_set):
        """promoted 和 pinned 的 h5ad 不出现在清理候选中。"""
        candidates = enumerate_cleanup_candidates(record_set)
        run_ids = {c.run_id for c in candidates}
        assert run_ids == {"superseded-1", "failed-1"}, (
            f"清理候选应仅含 superseded 和 failed run，实际: {run_ids}"
        )

    def test_recursive_scan_finds_nested_h5ad(self, record_set):
        candidates = enumerate_cleanup_candidates(record_set)
        superseded_paths = [c.path for c in candidates if c.run_id == "superseded-1"]
        path_names = {p.name for p in superseded_paths}
        assert path_names == {"output.h5ad", "deep.h5ad"}, (
            f"递归扫描应找到嵌套 subdir 中的 h5ad，实际: {path_names}"
        )

    def test_suffix_filter_excludes_txt(self, record_set):
        candidates = enumerate_cleanup_candidates(record_set)
        for c in candidates:
            assert c.path.suffix in LARGE_MATRIX_SUFFIXES, (
                f"清理候选不应包含非 h5ad 后缀文件: {c.path}"
            )

    def test_min_bytes_threshold(self, record_set):
        """min_bytes=1500：deep.h5ad(2048) 保留，output.h5ad(1024) 排除。"""
        candidates = enumerate_cleanup_candidates(record_set, min_bytes=1500)
        superseded = [c for c in candidates if c.run_id == "superseded-1"]
        assert len(superseded) == 1
        assert superseded[0].path.name == "deep.h5ad"
        assert superseded[0].size_bytes == 2048

    def test_files_not_deleted_after_enumeration(self, record_set):
        """铁律守卫：调用 enumerate_cleanup_candidates 后所有文件仍然存在。"""
        # 先收集所有候选路径
        candidates = enumerate_cleanup_candidates(record_set)
        for c in candidates:
            assert c.path.exists(), (
                f"enumerate_cleanup_candidates 禁止删除文件，但 {c.path} 不存在了"
            )

    def test_custom_suffixes(self, record_set):
        """自定义后缀：传入 (".txt",)，应只收集 .txt 文件。"""
        candidates = enumerate_cleanup_candidates(record_set, suffixes=(".txt",))
        assert len(candidates) >= 1
        for c in candidates:
            assert c.path.suffix == ".txt"

    def test_non_empty_proof_promoted_excluded(self, record_set):
        """非空跑证明：若类别过滤失效，promoted 的 h5ad 会混入结果。"""
        candidates = enumerate_cleanup_candidates(record_set)
        promoted_ids = {c.run_id for c in candidates if c.category == PROMOTED}
        assert len(promoted_ids) == 0, (
            f"promoted run 的 h5ad 不应出现在清理候选中，实际: {promoted_ids}"
        )
