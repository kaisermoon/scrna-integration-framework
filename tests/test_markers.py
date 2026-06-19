"""Tests for load_markers() -- synthetic CSV, no real data files."""

import os

import pytest

from scrna_integration.markers import load_markers

_SAMPLE = (
    "tissue,cell_type,marker,role,reference,notes\n"
    "tissue_A,cell_type_1,GENE_A,canonical,Author YEAR,core\n"
    "tissue_A,cell_type_1,GENE_B,canonical,Author YEAR,\n"
    "tissue_A,cell_type_1,GENE_E,optional,Author YEAR,\n"
    "tissue_A,cell_type_2,GENE_C,canonical,Author YEAR,\n"
    "tissue_A,cell_type_2,GENE_D,canonical,Author YEAR,\n"
    "tissue_A,cell_type_2,GENE_F,negative,Author YEAR,exclude\n"
)


@pytest.fixture
def markers_csv(tmp_path):
    p = tmp_path / "markers.csv"
    p.write_text(_SAMPLE)
    return str(p)


class TestLoadMarkers:
    def test_default_canonical_optional(self, markers_csv):
        path = markers_csv
        r = load_markers(path)
        assert "cell_type_1" in r
        assert "cell_type_2" in r
        assert "GENE_A" in r["cell_type_1"]
        assert "GENE_B" in r["cell_type_1"]
        assert "GENE_E" in r["cell_type_1"]
        assert "GENE_C" in r["cell_type_2"]
        assert "GENE_D" in r["cell_type_2"]
        # negative markers excluded by default
        assert "GENE_F" not in r["cell_type_2"]

    def test_negative_only(self, markers_csv):
        path = markers_csv
        r = load_markers(path, roles=("negative",))
        assert r["cell_type_2"] == ["GENE_F"]
        assert "cell_type_1" not in r  # no negative markers for cell_type_1

    def test_roles_none_three_layer(self, markers_csv):
        path = markers_csv
        r = load_markers(path, roles=None)
        assert r["cell_type_1"]["canonical"] == ["GENE_A", "GENE_B"]
        assert r["cell_type_1"]["optional"] == ["GENE_E"]
        assert r["cell_type_1"]["negative"] == []
        assert r["cell_type_2"]["canonical"] == ["GENE_C", "GENE_D"]
        assert r["cell_type_2"]["negative"] == ["GENE_F"]
        assert r["cell_type_2"]["optional"] == []

    def test_accepts_absolute_path(self, markers_csv):
        path = markers_csv
        r = load_markers(os.path.abspath(path))
        assert "cell_type_1" in r


_EMPTY_TEMPLATE = "tissue,cell_type,marker,role,reference,notes\n"


class TestLoadMarkersEmptyTemplate:
    """load_markers() on header-only template CSVs must not crash.

    Verifies that PI can safely load empty template CSVs (header-only,
    no data rows) without errors. This is the expected starting state
    before PI fills in real marker genes.
    """

    @pytest.fixture
    def empty_csv(self, tmp_path):
        p = tmp_path / "empty_template.csv"
        p.write_text(_EMPTY_TEMPLATE)
        return str(p)

    def test_empty_template_flat_returns_empty_dict(self, empty_csv):
        """Header-only template -> empty dict (no data rows to group)."""
        r = load_markers(empty_csv)
        assert r == {}

    def test_empty_template_roles_none_returns_empty_dict(self, empty_csv):
        """Three-layer mode on empty template -> empty dict."""
        r = load_markers(empty_csv, roles=None)
        assert r == {}

    def test_empty_template_negative_returns_empty_dict(self, empty_csv):
        """Negative-only filter on empty template -> empty dict."""
        r = load_markers(empty_csv, roles=("negative",))
        assert r == {}


# ---------------------------------------------------------------------------
# P1 xfail: load_markers roles=str（纯字符串而非 tuple）跨 pandas 版本脆弱
# ---------------------------------------------------------------------------


class TestLoadMarkersRolesBug:
    """确认 bug：roles 传纯字符串时缺乏类型守卫。"""

    @pytest.mark.xfail(
        reason=(
            "BUG: load_markers roles=str 时缺少类型守卫（markers.py:82）。"
            "pandas<2.0：isin(str) 按字符迭代静默返回 {}；"
            "pandas>=2.0：抛 TypeError 但报英文内部信息。"
            "修复方向：if isinstance(roles, str): raise TypeError('roles 必须为 list，不能传字符串')"
        ),
        strict=True,
    )
    def test_roles_str_raises_typeerror(self, markers_csv):
        """roles='canonical' 纯字符串应触发中文 TypeError，非依赖 pandas 内部报错。"""
        # 期望：产品代码加守卫后抛出中文 TypeError
        # 当前（无守卫）：pandas>=2.0 抛英文 TypeError（不匹配中文），
        #                pandas<2.0 静默返回 {}（不抛异常）。
        # 两种情况均不满足 match="必须为 list"，测试保持红色直到 bug 修复。
        with pytest.raises(TypeError, match="必须为 list"):
            load_markers(markers_csv, roles="canonical")
