"""PDB -> pdb2pqr (CHARMM) -> ESES surface -> CavityDetection -> PocketResult."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

from .result import (_res_sort, Pocket, PocketResult, read_obj, read_pqr, read_result_txt, voxel_from_log, write_xyzr)
from .structure import PreparedStructure, ResidueSpec, prepare_structure
from .tools import Tools, find_tools

ESES_SCRATCH = ("grid_info.txt", "intersection_info.txt", "partition_area.txt", "point_cloud.txt")


class PipelineError(RuntimeError):
    def __init__(self, step: str, message: str, log: str = ""):
        tail = "\n".join(log.strip().splitlines()[-25:])
        super().__init__(f"{step} failed: {message}" + (f"\n--- last lines of output ---\n{tail}" if tail else ""))
        self.step, self.log = step, log


def _run(cmd: List[str], cwd: str, step: str, timeout: Optional[float], verbose: bool,
         ok_codes: Sequence[int] = (0,)) -> str:
    if verbose:
        print(f"[ppdpy] {step}: {' '.join(cmd)}", file=sys.stderr, flush=True)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        raise PipelineError(step, f"timed out after {timeout} s", out)
    out = p.stdout or ""
    if p.returncode not in ok_codes:
        sig = ""
        if p.returncode < 0 or p.returncode > 128:
            n = -p.returncode if p.returncode < 0 else p.returncode - 128
            sig = f" (signal {n}{', illegal instruction: rebuild with the RemovePocket fix' if n == 4 else ''})"
        raise PipelineError(step, f"exit code {p.returncode}{sig}", out)
    if verbose:
        print(f"[ppdpy] {step}: done in {time.time() - t0:.1f} s", file=sys.stderr, flush=True)
    return f"$ {' '.join(cmd)}\n{out}\n"


# ------------------------------------------------------------------ single steps
def run_pdb2pqr(pdb_path: str, pqr_path: str, tools: Optional[Tools] = None, forcefield: str = "CHARMM",
                timeout: Optional[float] = None, verbose: bool = False) -> str:
    """pdb2pqr --ff=CHARMM --drop-water in.pdb out.pqr  (returns the program output)."""
    tools = tools or find_tools()
    cwd = os.path.dirname(os.path.abspath(pqr_path))
    log = _run(tools.pdb2pqr + [f"--ff={forcefield}", "--drop-water", os.path.abspath(pdb_path),
                                os.path.abspath(pqr_path)], cwd, "pdb2pqr", timeout, verbose)
    if not os.path.exists(pqr_path) or os.path.getsize(pqr_path) == 0:
        raise PipelineError("pdb2pqr", "no PQR written", log)
    return log


def run_eses(pqr_path: str, probe: float = 1.4, resolution: float = 0.4, margin: float = 1.0,
             tools: Optional[Tools] = None, timeout: Optional[float] = None, verbose: bool = False) -> Tuple[str, str]:
    """ESES MS_Intersection; writes marching_cubes.obj next to the PQR.  Returns (obj_path, output)."""
    tools = tools or find_tools(need_pdb2pqr=False)
    cwd = os.path.dirname(os.path.abspath(pqr_path))
    obj = os.path.join(cwd, "marching_cubes.obj")
    if os.path.exists(obj):
        os.remove(obj)
    log = _run(tools.eses + [os.path.basename(pqr_path), f"{probe:g}", f"{resolution:g}", f"{margin:g}"],
               cwd, "ESES", timeout, verbose)
    for f in ESES_SCRATCH:      # large intermediates, not needed afterwards
        try:
            os.remove(os.path.join(cwd, f))
        except OSError:
            pass
    if not os.path.exists(obj) or os.path.getsize(obj) == 0:
        raise PipelineError("ESES", "no marching_cubes.obj written", log)
    return obj, log


def _write_obj(verts, faces, path):
    with open(path, "w") as out:
        for v in np.asarray(verts, float):
            out.write(f"v {v[0]:g} {v[1]:g} {v[2]:g}\n")
        for f in np.asarray(faces, int):
            out.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


def run_cavity_detection(surface: Union[str, Tuple[np.ndarray, np.ndarray]],
                         atoms: Union[str, Tuple[np.ndarray, np.ndarray]],
                         min_depth: float = 4, min_area: float = 80, min_volume: float = 150,
                         tools: Optional[Tools] = None, workdir: Optional[str] = None,
                         timeout: Optional[float] = None, verbose: bool = False) -> Tuple[List[dict], Optional[float], str]:
    """Low level: run CavityDetection on a surface mesh and atom set.

    surface  path to an .obj (e.g. ESES marching_cubes.obj) or (vertices (V,3), faces (F,3) 0-based)
    atoms    path to an .xyzr file or (xyz (N,3), radii (N,))
    Returns (pockets, voxel_size_A, output).  Each pocket: dict(index, area, volume, depth_steps,
    depth (A), atoms = 0-based indices into the atom list)."""
    tools = tools or find_tools(need_pdb2pqr=False)
    own = workdir is None
    wd = tempfile.mkdtemp(prefix="ppd_cd_") if own else workdir
    try:
        obj = surface if isinstance(surface, str) else os.path.join(wd, "surface.obj")
        if not isinstance(surface, str):
            _write_obj(*surface, obj)
        xyzr = atoms if isinstance(atoms, str) else os.path.join(wd, "atoms.xyzr")
        if not isinstance(atoms, str):
            write_xyzr(atoms[0], atoms[1], xyzr)
        res = os.path.join(wd, "Result.txt")
        if os.path.exists(res):
            os.remove(res)
        # CavityDetection exits with status 1 on success
        log = _run(tools.cavity + [os.path.abspath(obj), os.path.abspath(xyzr), f"{min_depth:g}", f"{min_area:g}",
                                   f"{min_volume:g}"], wd, "CavityDetection", timeout, verbose, ok_codes=(0, 1))
        if not os.path.exists(res):
            raise PipelineError("CavityDetection", "no Result.txt written", log)
        pockets = read_result_txt(res)
        voxel = voxel_from_log(log)
        for p in pockets:
            p["depth"] = p["depth_steps"] * 0.5 * voxel if voxel else None
        return pockets, voxel, log
    finally:
        if own:
            shutil.rmtree(wd, ignore_errors=True)


# ------------------------------------------------------------------ full pipeline
def detect_pockets(structure: Union[str, os.PathLike, PreparedStructure],
                   chains: Optional[Iterable[str]] = None, peptide: Optional[Iterable[str]] = None,
                   model: int = 1, residues: ResidueSpec = None, peptide_max_len: int = 30,
                   probe: float = 1.4, resolution: float = 0.4, margin: float = 1.0,
                   min_depth: float = 4, min_area: float = 80, min_volume: float = 150,
                   workdir: Optional[str] = None, keep_files: bool = False, load_mesh: bool = False,
                   timeout: Optional[float] = None, verbose: bool = False, tools: Optional[Tools] = None,
                   cavity_detection: Optional[str] = None, eses: Optional[str] = None,
                   pdb2pqr: Optional[str] = None) -> PocketResult:
    """Detect pockets on a protein.

    structure  .pdb path (or PDB text, or a PreparedStructure), or a .pqr path (skips preparation
               and pdb2pqr; atoms are used as they are)
    chains     receptor chains (default: chains longer than peptide_max_len residues)
    peptide    chains treated as the bound peptide (removed; kept for PocketResult.evaluate())
    model      model number for NMR files;  residues: e.g. "20-80" to keep a receptor segment
    probe, resolution, margin           ESES surface settings (A)
    min_depth, min_area, min_volume     CavityDetection thresholds (steps, A^2, A^3)
    workdir    run in this folder (default: a temporary folder);  keep_files keeps it
    load_mesh  also return the ESES surface as PocketResult.mesh = (vertices, faces)
    """
    tools = tools or find_tools(cavity_detection, eses, pdb2pqr)
    settings = dict(probe=probe, resolution=resolution, margin=margin, min_depth=min_depth, min_area=min_area,
                    min_volume=min_volume, model=model, residues=residues if isinstance(residues, (str, type(None)))
                    else str(residues))
    own = workdir is None
    wd = tempfile.mkdtemp(prefix="ppd_") if own else os.path.abspath(workdir)
    os.makedirs(wd, exist_ok=True)
    log = ""
    t0 = time.time()
    try:
        src = os.fspath(structure) if isinstance(structure, os.PathLike) else structure
        is_pqr = isinstance(src, str) and "\n" not in src and src.lower().endswith(".pqr")
        pqr = os.path.join(wd, "receptor.pqr")
        if is_pqr:
            settings["input"] = os.path.abspath(src)
            shutil.copyfile(src, pqr)
            prep = None
        else:
            if isinstance(src, PreparedStructure):
                prep = src
            else:
                if isinstance(src, str) and "\n" not in src:
                    if not os.path.exists(src):
                        raise FileNotFoundError(src)
                    settings["input"] = os.path.abspath(src)
                prep = prepare_structure(src, chains=chains, peptide=peptide, model=model, residues=residues,
                                         peptide_max_len=peptide_max_len)
            rec_pdb = os.path.join(wd, "receptor.pdb")
            with open(rec_pdb, "w") as fh:
                fh.write(prep.receptor_pdb)
            if prep.peptide_pdb:
                with open(os.path.join(wd, "peptide.pdb"), "w") as fh:
                    fh.write(prep.peptide_pdb)
            log += run_pdb2pqr(rec_pdb, pqr, tools, timeout=timeout, verbose=verbose)

        atoms = read_pqr(pqr, prep.residue_order if prep else None)
        write_xyzr(atoms["xyz"], atoms["radius"], os.path.join(wd, "receptor.xyzr"))
        obj, out = run_eses(pqr, probe, resolution, margin, tools, timeout=timeout, verbose=verbose)
        log += out
        raw, voxel, out = run_cavity_detection(obj, os.path.join(wd, "receptor.xyzr"), min_depth, min_area,
                                               min_volume, tools, workdir=wd, timeout=timeout, verbose=verbose)
        log += out

        n = len(atoms["xyz"])
        pockets = []
        for p in raw:
            idx = np.array([i for i in p["atoms"] if 0 <= i < n], dtype=int)
            labels, seen = [], set()
            for i in idx:
                lab = f"{atoms['chain'][i]}:{atoms['resname'][i]}{atoms['resnum'][i]}{atoms['icode'][i]}"
                if lab not in seen:
                    seen.add(lab)
                    labels.append(lab)
            labels.sort(key=_res_sort)
            pockets.append(Pocket(index=p["index"], atom_indices=idx, area=p["area"], volume=p["volume"],
                                  depth_steps=p["depth_steps"], depth=p["depth"], residues=labels,
                                  centroid=atoms["xyz"][idx].mean(axis=0) if len(idx) else np.full(3, np.nan)))
        for r, p in enumerate(sorted(pockets, key=lambda q: -q.volume), 1):
            p.rank = r

        settings["runtime_s"] = round(time.time() - t0, 1)
        result = PocketResult(pockets=pockets, atoms=atoms, voxel_size=voxel, settings=settings,
                              receptor_chains=prep.receptor_chains if prep else sorted(set(atoms["chain"]) - {""}),
                              peptide_chains=prep.peptide_chains if prep else [],
                              peptide_xyz=prep.peptide_xyz if prep is not None and len(prep.peptide_xyz) else None,
                              notes=list(prep.notes) if prep else [], log=log,
                              workdir=wd if keep_files or not own else None,
                              mesh=read_obj(obj) if load_mesh else None)
        return result
    finally:
        if own and not keep_files:
            shutil.rmtree(wd, ignore_errors=True)
