"""Tests that do not need the compiled backend."""
from pathlib import Path

import numpy as np
import pytest

import ppdpy
from ppdpy.structure import parse_ranges

ROOT = Path(__file__).resolve().parents[1]
PDB = ROOT / "examples" / "1KNE.pdb"


def test_parse_ranges():
    assert parse_ranges("20-80") == [(20, 80)]
    assert parse_ranges("5-40, 45-80,90") == [(5, 40), (45, 80), (90, 90)]
    assert parse_ranges(range(20, 81)) == [(20, 80)]
    assert parse_ranges(None) == []


def test_prepare_splits_receptor_and_peptide():
    prep = ppdpy.prepare_structure(str(PDB))
    assert prep.receptor_chains == ["A"] and prep.peptide_chains == ["P"]
    assert len(prep.residue_order) == 52
    assert prep.peptide_xyz.shape == (50, 3)
    assert "HOH" not in prep.receptor_pdb
    assert any("HOH" in n for n in prep.notes)


def test_prepare_residue_range_and_errors():
    prep = ppdpy.prepare_structure(str(PDB), chains="A", residues="10-20")
    assert [int(r[1]) for r in prep.residue_order] == list(range(10, 21))
    with pytest.raises(ValueError):
        ppdpy.prepare_structure(str(PDB), chains="Z")


def test_read_result_txt():
    pockets = ppdpy.read_result_txt(str(ROOT / "tests" / "1KNE_Result.txt"))
    assert [p["index"] for p in pockets] == [0, 1]
    assert pockets[1]["volume"] == pytest.approx(424.421)
    assert pockets[1]["depth_steps"] == 10
    assert len(pockets[1]["atoms"]) == 87


def test_read_coordinates():
    xyz = ppdpy.read_coordinates(str(PDB), chains="P")
    assert xyz.shape == (50, 3) and np.isfinite(xyz).all()
