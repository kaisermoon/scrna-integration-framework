"""P1-a: SoupX R 脚本测试（独立文件，不改共享参数化测试）。

测试通过 subprocess 调用 Rscript，用最小合成夹具验证：
1. 成功路径产物 + status.json
2. tod 不裁 barcode（核心防回归）
3. barcode 重叠正常 + 完全不重叠 fail
4. 基因对齐 + 基因不一致
5. 估计失败 exit(1)（锁死"失败不冒充成功"红线）
6. 缺失输入 fail
7. 校正矩阵非负
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SOUPX_SCRIPT = str(ROOT / "scripts" / "soupx_run.R")

# ---------------------------------------------------------------------------
# 环境检测
# ---------------------------------------------------------------------------


def _r_available() -> tuple[str, bool, str]:
    """检测 Rscript + SoupX 是否可用。

    按 platform.py rscript_bin() 同逻辑：先从 CONDA_PREFIX 推导 scrna-integration-r
    环境路径，再回退 PATH 查找。

    Returns
    -------
    (rscript_path, available, reason)
    """
    candidates = []

    # 首选：从 CONDA_PREFIX 推导 R 环境中的 Rscript
    # （与 platform.py rscript_bin() 相同的推导逻辑）
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        envs_dir = os.path.dirname(conda_prefix)
        candidate = os.path.join(envs_dir, "scrna-integration-r", "bin", "Rscript")
        if os.path.isfile(candidate):
            candidates.append(candidate)

    # 回退：PATH 查找
    which_result = shutil.which("Rscript")
    if which_result is not None:
        candidates.append(which_result)
    # 始终保留 "Rscript" 作为最后尝试（可能通过 conda run 的环境注入生效）
    candidates.append("Rscript")

    for rscript in candidates:
        try:
            result = subprocess.run(
                [rscript, "--vanilla", "-e", "library(SoupX); library(Matrix)"],
                capture_output=True, text=True, timeout=30,
            )
            if result.returncode == 0:
                return (rscript, True, f"Rscript + SoupX available via {rscript}")
        except (subprocess.TimeoutExpired, OSError):
            continue

    return ("Rscript", False, "Rscript + SoupX not found (checked CONDA_PREFIX derivation + PATH)")


_RSCRIPT, _R_OK, _R_REASON = _r_available()
pytestmark = pytest.mark.skipif(not _R_OK, reason=_R_REASON)


# ---------------------------------------------------------------------------
# 夹具构造
# ---------------------------------------------------------------------------

def _write_10x_dir(
    out_dir: str,
    counts: sp.spmatrix,
    barcodes: list[str],
    gene_names: list[str],
    compress: bool = False,
) -> None:
    """将矩阵 + barcodes + features 写入类 10x 目录。

    counts 应为 genes(行) × cells(列) 的 CSC/CSR 稀疏矩阵。
    features.tsv 三列：feature_id, feature_name, feature_type，
    与 Python 侧 notebook 导出格式一致。
    """
    os.makedirs(out_dir, exist_ok=True)
    scipy.io.mmwrite(os.path.join(out_dir, "matrix.mtx"), counts)
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        f.write("\n".join(barcodes) + "\n")
    with open(os.path.join(out_dir, "features.tsv"), "w") as f:
        for gn in gene_names:
            f.write(f"{gn}\t{gn}\tGene Expression\n")
    # compress 模式下 gzip，R 侧 read_10x_like 兼容 .gz
    if compress:
        import gzip as _gz
        for stem in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
            fp = os.path.join(out_dir, stem)
            if os.path.exists(fp):
                with open(fp, "rb") as src, _gz.open(fp + ".gz", "wb") as dst:
                    dst.writelines(src)
                os.remove(fp)


def _make_nice_fixture(
    tmp: Path,
    n_genes: int = 150,
    n_cells: int = 30,
    n_droplets: int = 80,
    rng: np.random.Generator | None = None,
) -> tuple[str, str, list[str], list[str]]:
    """构造正常夹具：三群细胞 + 低 UMI 空液滴，SoupX autoEstCont 可成功。

    关键约束（SoupX 1.6.2）:
    1. estimateSoup soupRange=c(0,100)：空液滴 nDropUMIs < 100
    2. soupQuantile=0.9：候选 marker genes 需 soupProf$est 在 top 10%
       → marker genes 的空液滴表达须高于 background genes
       （反向设计：marker 基因在空液滴中背景更高，因为其在 soup 中丰度大）
    3. 统计学效力：需足够 cells (>=20) + genes (>=100) 使超几何检验显著

    基因分组:
    - GENE 1-20:   Type A markers (cells 0-9)
    - GENE 21-40:  Type B markers (cells 10-19)
    - GENE 41-60:  Type C markers (cells 20-29)
    - GENE 61-150: background genes (all cells express uniformly low)

    Returns
    -------
    (filtered_dir, raw_dir, cell_barcodes, gene_names)
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n_per_type = n_cells // 3
    n_type_a = n_per_type
    n_type_b = n_per_type

    n_marker_per_type = 20
    n_marker_a = n_marker_per_type
    n_marker_b = n_marker_per_type
    n_marker_c = n_marker_per_type

    gene_names = [f"GENE{i:04d}" for i in range(1, n_genes + 1)]
    cell_barcodes = [f"CELL{i:06d}-1" for i in range(1, n_cells + 1)]
    droplet_barcodes = list(cell_barcodes)
    for i in range(n_cells, n_droplets):
        droplet_barcodes.append(f"DROPLET{i:06d}-1")

    # Gene index ranges
    idx_a = (0, n_marker_a)
    idx_b = (n_marker_a, n_marker_a + n_marker_b)
    idx_c = (n_marker_a + n_marker_b, n_marker_a + n_marker_b + n_marker_c)
    idx_bg = n_marker_a + n_marker_b + n_marker_c
    _all_markers = [(0, n_type_a, idx_a), (n_type_a, n_type_a + n_type_b, idx_b),
                    (n_type_a + n_type_b, n_cells, idx_c)]

    # Build filtered matrix: genes × cells
    filt_data = np.zeros((n_genes, n_cells), dtype=np.float64)
    for c in range(n_cells):
        for _cstart, _cend, _idx in _all_markers:
            if _cstart <= c < _cend:
                ms, me = _idx
                break
        else:
            ms, me = 0, 0
        for g in range(n_genes):
            if ms <= g < me:
                filt_data[g, c] = max(1, int(rng.poisson(rng.uniform(80, 200))))
            elif g < idx_bg:
                filt_data[g, c] = 0  # other clusters' markers off
            else:
                filt_data[g, c] = max(1, int(rng.poisson(rng.uniform(0.5, 2))))

    # Build raw droplet matrix: genes × droplets
    raw_data = np.zeros((n_genes, n_droplets), dtype=np.float64)
    for d in range(n_droplets):
        if d < n_cells:
            # Cell droplet: same structure as filtered
            for _cstart, _cend, _idx in _all_markers:
                if _cstart <= d < _cend:
                    ms, me = _idx
                    break
            else:
                ms, me = 0, 0
            for g in range(n_genes):
                if ms <= g < me:
                    raw_data[g, d] = max(1, int(rng.poisson(rng.uniform(80, 200))))
                elif g < idx_bg:
                    raw_data[g, d] = 0
                else:
                    raw_data[g, d] = max(1, int(rng.poisson(rng.uniform(0.5, 2))))
        else:
            # Empty droplet: marker genes have higher soup presence than bg
            # (soupQuantile=0.9 选 soupProf top 10%，marker 需高于 bg)
            for g in range(idx_bg):
                raw_data[g, d] = rng.poisson(rng.uniform(0.3, 1.2))
            for g in range(idx_bg, n_genes):
                raw_data[g, d] = rng.poisson(rng.uniform(0.05, 0.3))

    filtered_counts = sp.csc_matrix(filt_data)
    raw_counts = sp.csc_matrix(raw_data)

    filtered_dir = str(tmp / "filtered")
    raw_dir = str(tmp / "raw")
    _write_10x_dir(filtered_dir, filtered_counts, cell_barcodes, gene_names)
    _write_10x_dir(raw_dir, raw_counts, droplet_barcodes, gene_names)

    return filtered_dir, raw_dir, cell_barcodes, gene_names


def _make_degenerate_fixture(
    tmp: Path,
    n_genes: int = 50,
    n_cells: int = 1,
    n_droplets: int = 20,
) -> tuple[str, str, list[str], list[str]]:
    """构造退化夹具：1 个全零 cell + 低 UMI 空液滴，autoEstCont 必失败。

    单 cell 全零 filtered + 空液滴 UMI < 50（确保 soupRange 命中，不 NaN）。
    autoEstCont 因无法找到 marker genes 而抛错。
    """
    gene_names = [f"GENE{i:04d}" for i in range(1, n_genes + 1)]
    cell_barcodes = ["CELL000001-1"]
    droplet_barcodes = ["CELL000001-1"] + [
        f"DROPLET{i:06d}-1" for i in range(1, n_droplets)
    ]

    # filtered: single cell with all zeros — 无表达，autoEstCont 必失败
    filt_data = np.zeros((n_genes, n_cells), dtype=np.float64)
    filtered_counts = sp.csc_matrix(filt_data)

    # raw: empty droplets have very low UMI (< 50)，cell droplet also low
    rng = np.random.default_rng(99)
    raw_data = np.zeros((n_genes, n_droplets), dtype=np.float64)
    for d in range(n_droplets):
        for g in range(n_genes):
            raw_data[g, d] = rng.poisson(0.5 if d > 0 else 0.3)
    # Ensure cell droplet also < 100 UMI
    raw_counts = sp.csc_matrix(raw_data)

    raw_dir = str(tmp / "raw_degen")
    filtered_dir = str(tmp / "filtered_degen")
    _write_10x_dir(filtered_dir, filtered_counts, cell_barcodes, gene_names)
    _write_10x_dir(raw_dir, raw_counts, droplet_barcodes, gene_names)

    return filtered_dir, raw_dir, cell_barcodes, gene_names


# ---------------------------------------------------------------------------
# 辅助函数
# ---------------------------------------------------------------------------

def _run_soupx(
    work_dir: str,
    filtered_dir: str,
    raw_dir: str,
    sample_id: str = "test",
    timeout: int = 120,
) -> subprocess.CompletedProcess:
    """调用 Rscript soupx_run.R 并返回 CompletedProcess。"""
    return subprocess.run(
        [_RSCRIPT, "--vanilla", SOUPX_SCRIPT,
         os.path.abspath(work_dir), os.path.abspath(filtered_dir),
         os.path.abspath(raw_dir), sample_id],
        capture_output=True, text=True, timeout=timeout,
    )


def _read_status(work_dir: str) -> dict:
    """读 soupx_status.json。"""
    with open(os.path.join(work_dir, "soupx_status.json")) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# 测试用例
# ---------------------------------------------------------------------------

class TestSoupXSuccess:
    """正常夹具：成功路径验证。"""

    def test_success_writes_corrected_and_status(self, tmp_path: Path) -> None:
        """正常夹具 → exit 0；产物齐全 + status=success。"""
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)

        assert result.returncode == 0, f"exit={result.returncode}, stderr={result.stderr[:300]}"

        # corrected_counts.mtx 存在且 shape == genes × cells
        mtx_path = os.path.join(work_dir, "corrected_counts.mtx")
        assert os.path.exists(mtx_path), "corrected_counts.mtx 未产出"
        corrected = sp.csr_matrix(scipy.io.mmread(mtx_path))
        n_genes_expected = len(gene_names)
        n_cells_expected = len(cell_barcodes)
        assert corrected.shape == (n_genes_expected, n_cells_expected), \
            f"corrected shape={corrected.shape}, expected=({n_genes_expected},{n_cells_expected})"

        # barcodes.tsv 行序 == 输入 filtered barcode
        bc_path = os.path.join(work_dir, "barcodes.tsv")
        assert os.path.exists(bc_path)
        with open(bc_path) as f:
            out_barcodes = [line.strip() for line in f if line.strip()]
        assert out_barcodes == cell_barcodes, "barcode 顺序不匹配"

        # features.tsv 第 2 列 == 输入基因序
        feat_path = os.path.join(work_dir, "features.tsv")
        assert os.path.exists(feat_path)
        with open(feat_path) as f:
            out_features = [line.strip().split("\t") for line in f if line.strip()]
        out_gene_names = [cols[1] for cols in out_features]
        assert out_gene_names == gene_names, "基因顺序不匹配"

        # soupx_status.json status=success
        status = _read_status(work_dir)
        assert status["status"] == "success"
        assert 0.0 <= status["rho_global"] < 1.0
        assert status["n_cells_filtered"] == n_cells_expected
        assert status["n_genes"] == n_genes_expected
        assert status["n_common_genes"] == n_genes_expected
        assert status["method_params"] in ("default", "relaxed", "forceAccept")

        # rho.txt 存在
        rho_path = os.path.join(work_dir, "rho.txt")
        assert os.path.exists(rho_path)

    def test_tod_not_trimmed_to_filtered(self, tmp_path: Path) -> None:
        """断言 R 侧日志记录的 droplet 数 == 完整 raw D，而非 N。

        核心防回归：tod 未被裁到 filtered barcode 交集。
        """
        n_cells, n_droplets = 30, 80
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(
            tmp_path, n_cells=n_cells, n_droplets=n_droplets)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        status = _read_status(work_dir)
        # status.json 记录 n_droplets == 完整 raw droplet 数（80），
        # 而非 filtered 数（30）
        assert status["n_droplets"] == n_droplets, \
            f"n_droplets={status['n_droplets']}，应为 {n_droplets}（完整 raw）；" \
            f"若为 {n_cells} 则说明 tod 被裁"

        # 同时日志中应出现 {n_droplets} droplets 的计数
        stdout = result.stdout
        assert "tod=" in stdout, f"stdout 缺 tod= 日志: {stdout[:500]}"
        # 提取 tod 日志中的 droplet 数
        for line in stdout.split("\n"):
            if "tod=" in line and "droplets" in line:
                assert f"{n_droplets} droplets" in line or f"x {n_droplets}" in line, \
                    f"tod droplet count mismatch in log: {line}"

    def test_barcode_superset_not_intersected(self, tmp_path: Path) -> None:
        """filtered barcode 全部 in raw，结果 cells 数 == N，不因交集裁剪丢细胞。"""
        n_cells = 30
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(
            tmp_path, n_cells=n_cells)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        status = _read_status(work_dir)
        assert status["n_cells_filtered"] == n_cells, \
            f"cells={status['n_cells_filtered']}，应为 {n_cells}（未因交集裁剪丢细胞）"

    def test_corrected_counts_nonneg(self, tmp_path: Path) -> None:
        """成功路径校正矩阵所有值 >=0（SoupX 输出应为非负整数计数）。"""
        filtered_dir, raw_dir, _, _ = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        mtx_path = os.path.join(work_dir, "corrected_counts.mtx")
        corrected = scipy.io.mmread(mtx_path)
        # 检查所有非零值非负
        if sp.issparse(corrected):
            assert np.all(corrected.data >= 0), "校正矩阵含负值"
        else:
            assert np.all(corrected >= 0), "校正矩阵含负值"


class TestSoupXGeneAlignment:
    """基因对齐边界情况。"""

    def test_gene_mismatch_aligns_genes_only(self, tmp_path: Path) -> None:
        """raw 比 filtered 多若干基因 → 只按 filtered 基因集对齐，cells 不变。"""
        n_cells = 21
        n_droplets_empty = 40
        n_droplets = n_cells + n_droplets_empty

        rng = np.random.default_rng(123)
        # filtered 基因集：100 genes (GENE001-GENE100)
        # raw 基因集：150 genes (GENE051-GENE200)，交集 = GENE051-GENE100 = 50 genes
        gene_names_filt = [f"GENE{i:03d}" for i in range(1, 101)]
        gene_names_raw = [f"GENE{i:03d}" for i in range(51, 201)]
        common_names = sorted(set(gene_names_filt) & set(gene_names_raw))
        n_common = len(common_names)

        cell_barcodes = [f"CELL{i:06d}-1" for i in range(1, n_cells + 1)]
        droplet_barcodes = list(cell_barcodes) + [
            f"DROPLET{i:06d}-1" for i in range(1, n_droplets_empty + 1)
        ]

        # filtered: structured 3-group expression（需足够细胞让 autoEstCont 成功）
        filt_data = np.zeros((len(gene_names_filt), n_cells), dtype=np.float64)
        n_per = n_cells // 3
        for c in range(n_cells):
            if c < n_per:
                _mrk = list(range(0, 20))
            elif c < 2 * n_per:
                _mrk = list(range(20, 40))
            else:
                _mrk = list(range(40, 60))
            for g in range(len(gene_names_filt)):
                if g in _mrk:
                    filt_data[g, c] = max(1, int(rng.poisson(rng.uniform(80, 200))))
                elif g < 60:
                    filt_data[g, c] = 0
                else:
                    filt_data[g, c] = max(1, int(rng.poisson(rng.uniform(0.5, 2))))

        # raw: cell droplets similar, empty droplets with markers higher than bg
        raw_data = np.zeros((len(gene_names_raw), n_droplets), dtype=np.float64)
        for d in range(n_droplets):
            if d < n_cells:
                if d < n_per:
                    _mrk = [gene_names_raw.index(n) for n in common_names if n in gene_names_raw[:20]]
                elif d < 2 * n_per:
                    _mrk = [gene_names_raw.index(n) for n in common_names if n in gene_names_raw[20:40]]
                else:
                    _mrk = [gene_names_raw.index(n) for n in common_names if n in gene_names_raw[40:60]]
                _mrk_set = set(_mrk)
                for g in range(len(gene_names_raw)):
                    if g in _mrk_set:
                        raw_data[g, d] = max(1, int(rng.poisson(rng.uniform(80, 200))))
                    elif g < 60:
                        raw_data[g, d] = 0
                    else:
                        raw_data[g, d] = max(1, int(rng.poisson(rng.uniform(0.5, 2))))
            else:
                # Empty droplet: markers in top 60 genes have higher bg
                for g in range(min(60, len(gene_names_raw))):
                    raw_data[g, d] = rng.poisson(rng.uniform(0.3, 1.2))
                for g in range(60, len(gene_names_raw)):
                    raw_data[g, d] = rng.poisson(rng.uniform(0.05, 0.3))

        filtered_dir = str(tmp_path / "filtered_mismatch")
        raw_dir = str(tmp_path / "raw_mismatch")
        _write_10x_dir(filtered_dir, sp.csc_matrix(filt_data), cell_barcodes, gene_names_filt)
        _write_10x_dir(raw_dir, sp.csc_matrix(raw_data), droplet_barcodes, gene_names_raw)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0, \
            f"exit={result.returncode}, stderr={result.stderr[:200]}"

        status = _read_status(work_dir)
        assert status["n_common_genes"] == n_common, \
            f"common_genes={status['n_common_genes']}, expected={n_common}"
        assert status["n_cells_filtered"] == n_cells

        # corrected_counts 基因数 = common_genes，cells = n_cells
        mtx_path = os.path.join(work_dir, "corrected_counts.mtx")
        corrected = sp.csr_matrix(scipy.io.mmread(mtx_path))
        assert corrected.shape == (n_common, n_cells), \
            f"corrected shape={corrected.shape}, expected=({n_common},{n_cells})"

        # barcodes 顺序不变
        with open(os.path.join(work_dir, "barcodes.tsv")) as f:
            out_bc = [line.strip() for line in f if line.strip()]
        assert out_bc == cell_barcodes

        # features 第 2 列为共同基因（按 filtered 基因序排列）
        with open(os.path.join(work_dir, "features.tsv")) as f:
            out_features = [line.strip().split("\t") for line in f if line.strip()]
        out_genes = [cols[1] for cols in out_features]
        expected_common = [g for g in gene_names_filt if g in set(gene_names_raw)]
        assert out_genes == expected_common, \
            f"output gene order mismatch (len={len(out_genes)} vs {len(expected_common)})"

    def test_gene_mismatch_common_less_than_10_fails(self, tmp_path: Path) -> None:
        """共同基因 < 10 → exit 1 + status=failed。"""
        gene_names_filt = [f"GENE_A{i:04d}" for i in range(1, 6)]
        gene_names_raw = [f"GENE_B{i:04d}" for i in range(1, 6)]
        cell_barcodes = ["CELL000001-1"]
        droplet_barcodes = ["CELL000001-1", "DROPLET000001-1"]

        rng = np.random.default_rng(42)
        filt_data = rng.poisson(5, (5, 1)).astype(np.float64)
        raw_data = rng.poisson(3, (5, 2)).astype(np.float64)

        filtered_dir = str(tmp_path / "filtered_few")
        raw_dir = str(tmp_path / "raw_few")
        _write_10x_dir(filtered_dir, sp.csc_matrix(filt_data), cell_barcodes, gene_names_filt)
        _write_10x_dir(raw_dir, sp.csc_matrix(raw_data), droplet_barcodes, gene_names_raw)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 1, f"expected exit 1, got {result.returncode}"

        status = _read_status(work_dir)
        assert status["status"] == "failed"
        assert "common" in status["reason"].lower() or "few" in status["reason"].lower()
        assert not os.path.exists(os.path.join(work_dir, "corrected_counts.mtx"))


class TestSoupXFailure:
    """失败路径硬退出验证。"""

    def test_no_overlap_barcodes_fails(self, tmp_path: Path) -> None:
        """filtered barcode 与 raw 完全不相交 → exit 1 + status=failed。"""
        gene_names = [f"GENE{i:04d}" for i in range(1, 51)]
        cell_barcodes = ["CELL000001-1", "CELL000002-1"]
        droplet_barcodes = ["DROPLET000001-1", "DROPLET000002-1"]

        rng = np.random.default_rng(1)
        filt_data = rng.poisson(10, (50, 2)).astype(np.float64)
        raw_data = rng.poisson(5, (50, 2)).astype(np.float64)

        filtered_dir = str(tmp_path / "filtered_no_overlap")
        raw_dir = str(tmp_path / "raw_no_overlap")
        _write_10x_dir(filtered_dir, sp.csc_matrix(filt_data), cell_barcodes, gene_names)
        _write_10x_dir(raw_dir, sp.csc_matrix(raw_data), droplet_barcodes, gene_names)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 1, f"expected exit 1, got {result.returncode}"

        status = _read_status(work_dir)
        assert status["status"] == "failed"
        assert "overlap" in status["reason"].lower() or "barcode" in status["reason"].lower()
        # 绝不产出 corrected_counts.mtx
        assert not os.path.exists(os.path.join(work_dir, "corrected_counts.mtx")), \
            "失败路径不应产出 corrected_counts.mtx"

    def test_estimation_failure_exits_nonzero(self, tmp_path: Path) -> None:
        """退化夹具（单细胞全零）→ autoEstCont 失败 → exit 1 + 不产出 corrected。

        这条锁死"失败不冒充成功"红线。
        """
        filtered_dir, raw_dir, _, _ = _make_degenerate_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir, timeout=120)
        assert result.returncode == 1, \
            f"expected exit 1 for estimation failure, got {result.returncode}"

        status = _read_status(work_dir)
        assert status["status"] == "failed", \
            f"expected status=failed, got {status['status']}"
        assert "autoEstCont" in status["reason"] or "fail" in status["reason"].lower(), \
            f"reason should mention estimation failure, got: {status['reason']}"

        # 红线：绝不产出 corrected_counts.mtx
        assert not os.path.exists(os.path.join(work_dir, "corrected_counts.mtx")), \
            "估计失败绝不应产出 corrected_counts.mtx（冒充成功红线）"

    def test_missing_input_dir_fails(self, tmp_path: Path) -> None:
        """raw 目录缺 matrix 文件 → exit 1 + stderr 明确原因。"""
        # filtered 正常，raw 目录不存在
        filtered_dir, _, cell_barcodes, gene_names = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")
        missing_raw = str(tmp_path / "nonexistent")

        result = _run_soupx(work_dir, filtered_dir, missing_raw)
        assert result.returncode == 1, f"expected exit 1, got {result.returncode}"

        status = _read_status(work_dir)
        assert status["status"] == "failed"
        assert not os.path.exists(os.path.join(work_dir, "corrected_counts.mtx"))

    def test_missing_filtered_input_fails(self, tmp_path: Path) -> None:
        """filtered 目录缺 matrix 文件 → exit 1。"""
        _, raw_dir, _, _ = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")
        missing_filtered = str(tmp_path / "nonexistent_filtered")

        result = _run_soupx(work_dir, missing_filtered, raw_dir)
        assert result.returncode == 1

        status = _read_status(work_dir)
        assert status["status"] == "failed"


class TestSoupXEdgeCases:
    """边界情况。"""

    def test_corrected_barcodes_match_filtered_order(self, tmp_path: Path) -> None:
        """明确验证 barcodes.tsv 行序严格 == 输入 filtered barcode。

        notebook 侧逐行比对，顺序不等即跳过整样本。
        """
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        bc_path = os.path.join(work_dir, "barcodes.tsv")
        with open(bc_path) as f:
            out_barcodes = [line.strip() for line in f if line.strip()]
        # 逐元素比对
        for i, (expected, actual) in enumerate(zip(cell_barcodes, out_barcodes, strict=True)):
            assert expected == actual, \
                f"barcode[{i}] mismatch: expected={expected}, actual={actual}"

    def test_features_col2_matches_gene_order(self, tmp_path: Path) -> None:
        """明确验证 features.tsv 第 2 列（iloc[:,1]）严格 == 输入基因序。"""
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        feat_path = os.path.join(work_dir, "features.tsv")
        with open(feat_path) as f:
            out_genes = [line.strip().split("\t")[1] for line in f if line.strip()]
        assert out_genes == gene_names, \
            "features.tsv col2 gene order mismatch"

    def test_rho_file_contains_valid_number(self, tmp_path: Path) -> None:
        """rho.txt 包含合法的 0-1 数字。"""
        filtered_dir, raw_dir, _, _ = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")

        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        rho_path = os.path.join(work_dir, "rho.txt")
        with open(rho_path) as f:
            rho_val = float(f.read().strip())
        assert 0.0 <= rho_val < 1.0, f"rho={rho_val} out of [0,1) range"

    def test_raw_compressed_input_works(self, tmp_path: Path) -> None:
        """.gz 压缩 raw input 正常读取。"""
        # 用正常夹具但 raw 侧加 .gz 包装
        filtered_dir, raw_dir, cell_barcodes, gene_names = _make_nice_fixture(tmp_path)
        # 对 raw 目录加 .gz
        for stem in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
            fp_src = os.path.join(raw_dir, stem)
            fp_dst = fp_src + ".gz"
            import gzip
            with open(fp_src, "rb") as src, gzip.open(fp_dst, "wb") as dst:
                dst.writelines(src)
            os.remove(fp_src)

        work_dir = str(tmp_path / "work_gz")
        result = _run_soupx(work_dir, filtered_dir, raw_dir)
        assert result.returncode == 0

        status = _read_status(work_dir)
        assert status["status"] == "success"

    def test_sample_id_in_stdout(self, tmp_path: Path) -> None:
        """sample_id 出现在日志中（用于日志标识）。"""
        filtered_dir, raw_dir, _, _ = _make_nice_fixture(tmp_path)
        work_dir = str(tmp_path / "work")
        sample_id = "my_sample_42"

        result = _run_soupx(work_dir, filtered_dir, raw_dir, sample_id=sample_id)
        assert result.returncode == 0
        assert sample_id in result.stdout, \
            f"sample_id '{sample_id}' not found in stdout"


def test_module_level_skip_when_no_r() -> None:
    """模块级 pytestmark 覆盖：R 不可用时所有测试 skip。"""
    # 此测试本身不依赖 R——如果 skip mark 生效，其他测试已被跳过。
    # 这里确保 _R_REASON 是合法字符串。
    assert isinstance(_R_REASON, str) and len(_R_REASON) > 0
