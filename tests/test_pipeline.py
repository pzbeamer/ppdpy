"""End-to-end test (about a minute). Skipped unless pdb2pqr, ESES and CavityDetection are found."""
from pathlib import Path

import pytest

import ppdpy

ROOT = Path(__file__).resolve().parents[1]

try:
    ppdpy.find_tools()
    HAVE_TOOLS = True
except ppdpy.ToolNotFound:
    HAVE_TOOLS = False

pytestmark = [pytest.mark.slow, pytest.mark.skipif(not HAVE_TOOLS, reason="backend not built (ppd build-backend)")]


def test_1kne_matches_reference(tmp_path):
    res = ppdpy.detect_pockets(str(ROOT / "examples" / "1KNE.pdb"))
    ref = {p["index"]: p for p in ppdpy.read_result_txt(str(ROOT / "tests" / "1KNE_Result.txt"))}
    assert len(res) == len(ref) == 2
    for p in res:
        r = ref[p.index]
        assert p.volume == pytest.approx(r["volume"], rel=1e-3)
        assert p.area == pytest.approx(r["area"], rel=1e-3)
        assert p.depth_steps == r["depth_steps"]
        assert p.atom_indices.tolist() == r["atoms"]
    best = res.evaluate()[0]
    assert best["pocket"] == 1
    assert best["residue_recall"] > 0.9 and best["dcc"] < 5
    res.save(str(tmp_path))
    assert (tmp_path / "pockets.json").exists()
