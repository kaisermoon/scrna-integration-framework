"""Tests for load_markers() — synthetic CSV, no real data files."""

import os

import pytest

from scrna_integration.markers import load_markers

_SAMPLE = (
    "tissue,cell_type,marker,role,reference,notes\n"
    "gastric,SPEM,TFF2,canonical,Goldenring 2017,core\n"
    "gastric,SPEM,MUC6,canonical,Goldenring 2017,\n"
    "gastric,SPEM,CD44,optional,Ref 2020,\n"
    "gastric,pit_cell,MUC5AC,canonical,Nowicki 2023,\n"
    "gastric,pit_cell,TFF1,canonical,Nowicki 2023,\n"
    "gastric,pit_cell,CD44v9,negative,Ref 2021,exclude\n"
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
        assert "SPEM" in r
        assert "pit_cell" in r
        assert "TFF2" in r["SPEM"]
        assert "MUC6" in r["SPEM"]
        assert "CD44" in r["SPEM"]
        assert "MUC5AC" in r["pit_cell"]
        assert "TFF1" in r["pit_cell"]
        # negative markers excluded by default
        assert "CD44v9" not in r["pit_cell"]

    def test_negative_only(self, markers_csv):
        path = markers_csv
        r = load_markers(path, roles=("negative",))
        assert r["pit_cell"] == ["CD44v9"]
        assert "SPEM" not in r  # no negative markers for SPEM

    def test_roles_none_three_layer(self, markers_csv):
        path = markers_csv
        r = load_markers(path, roles=None)
        assert r["SPEM"]["canonical"] == ["TFF2", "MUC6"]
        assert r["SPEM"]["optional"] == ["CD44"]
        assert r["SPEM"]["negative"] == []
        assert r["pit_cell"]["canonical"] == ["MUC5AC", "TFF1"]
        assert r["pit_cell"]["negative"] == ["CD44v9"]
        assert r["pit_cell"]["optional"] == []

    def test_accepts_absolute_path(self, markers_csv):
        path = markers_csv
        r = load_markers(os.path.abspath(path))
        assert "SPEM" in r
