"""ppdpy - Python wrapper around PPD (CavityDetection) + ESES + PDB2PQR.

PPD: Zhao R, Cang Z, Tong Y, Wei G-W. Protein pocket detection via convex hull surface
evolution and associated Reeb graph. Bioinformatics 34(17):i830-i837 (2018).
https://github.com/rdzhao/PPD  -- the pocket detection itself is entirely their work.

    import ppdpy
    res = ppdpy.detect_pockets("1KNE.pdb", chains="A")
    for p in res.by_volume():
        print(p)
    print(res.evaluate()[0])        # if the file contained a bound peptide
"""
from .pipeline import PipelineError, detect_pockets, run_cavity_detection, run_eses, run_pdb2pqr
from .result import Pocket, PocketResult, read_obj, read_pqr, read_result_txt
from .structure import PreparedStructure, prepare_structure, read_coordinates
from .tools import ToolNotFound, Tools, backend_home, doctor, find_tools

__version__ = "0.1.0"
__all__ = ["detect_pockets", "run_cavity_detection", "run_eses", "run_pdb2pqr", "PipelineError",
           "Pocket", "PocketResult", "read_obj", "read_pqr", "read_result_txt",
           "PreparedStructure", "prepare_structure", "read_coordinates",
           "ToolNotFound", "Tools", "backend_home", "doctor", "find_tools"]
