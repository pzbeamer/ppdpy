"""Locate the three external programs the pipeline needs.

Search order for each tool:
  1. the explicit path passed to find_tools() / detect_pockets()
  2. environment variable  PPD_CAVITYDETECTION / PPD_ESES / PPD_PDB2PQR
  3. the PATH  (CavityDetection, MS_Intersection, pdb2pqr30 / pdb2pqr)
  4. $PPDPY_HOME/bin (default ~/.local/share/ppdpy/bin), where `ppd build-backend` installs
  5. a few common manual build locations (see _CANDIDATES)
  6. for pdb2pqr only: "python -m pdb2pqr" if the module is importable
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


class ToolNotFound(RuntimeError):
    pass


def backend_home() -> Path:
    """Install location of the compiled backend ($PPDPY_HOME, default ~/.local/share/ppdpy)."""
    return Path(os.path.expanduser(os.environ.get("PPDPY_HOME", "~/.local/share/ppdpy")))


_CANDIDATES = {
    "cavity": [
        "~/PPD/build/CavityDetection", "~/PPD-master/build/CavityDetection", "~/PPD-master/CavityDetection",
        "~/PPD-bundle/PPD-master/build/CavityDetection",
    ],
    "eses": ["~/ESES/build/MS_Intersection", "~/PPD-bundle/ESES/build/MS_Intersection"],
    "pdb2pqr": [
        f"~/{base}/envs/{env}/bin/{exe}"
        for base in ("miniforge3", "miniconda3", "anaconda3", "mambaforge", "opt/anaconda3", "opt/miniconda3")
        for env in ("ppd", "pdb2pqr")
        for exe in ("pdb2pqr30", "pdb2pqr")
    ],
}
_ENV = {"cavity": "PPD_CAVITYDETECTION", "eses": "PPD_ESES", "pdb2pqr": "PPD_PDB2PQR"}
_NAMES = {"cavity": ["CavityDetection"], "eses": ["MS_Intersection"], "pdb2pqr": ["pdb2pqr30", "pdb2pqr"]}


def _is_exe(p) -> bool:
    return p is not None and os.path.isfile(p) and os.access(p, os.X_OK)


def _home_candidates(key):
    return [backend_home() / "bin" / n for n in _NAMES[key]] if key != "pdb2pqr" else []


def _find(key, explicit=None) -> List[str]:
    """Return the command (list of argv words) for tool `key`."""
    if explicit:
        p = os.path.expanduser(str(explicit))
        if _is_exe(p):
            return [p]
        w = shutil.which(p)
        if w:
            return [w]
        raise ToolNotFound(f"{key}: {explicit} is not an executable file")
    env = os.environ.get(_ENV[key])
    if env:
        p = os.path.expanduser(env)
        if _is_exe(p) or shutil.which(p):
            return [shutil.which(p) or p]
        raise ToolNotFound(f"{_ENV[key]}={env} is not an executable file")
    for name in _NAMES[key]:
        w = shutil.which(name)
        if w:
            return [w]
    for c in list(_home_candidates(key)) + [Path(os.path.expanduser(c)) for c in _CANDIDATES[key]]:
        if _is_exe(str(c)):
            return [str(Path(c).resolve())]
    if key == "pdb2pqr" and importlib.util.find_spec("pdb2pqr") is not None:
        return [sys.executable, "-m", "pdb2pqr"]
    hint = {
        "cavity": "run `ppd build-backend` (or set PPD_CAVITYDETECTION=/path/to/CavityDetection)",
        "eses": "run `ppd build-backend` (or set PPD_ESES=/path/to/MS_Intersection)",
        "pdb2pqr": "pip install pdb2pqr  (or set PPD_PDB2PQR=/path/to/pdb2pqr30)",
    }[key]
    raise ToolNotFound(f"could not find {'/'.join(_NAMES[key])}: {hint}")


@dataclass
class Tools:
    cavity: List[str]
    eses: List[str]
    pdb2pqr: List[str]

    def describe(self) -> str:
        return "\n".join(f"{k:8s} {' '.join(v)}" for k, v in
                         (("PPD", self.cavity), ("ESES", self.eses), ("pdb2pqr", self.pdb2pqr)))


def find_tools(cavity: Optional[str] = None, eses: Optional[str] = None,
               pdb2pqr: Optional[str] = None, need_pdb2pqr: bool = True) -> Tools:
    """Resolve the three programs (raises ToolNotFound with a hint if one is missing)."""
    return Tools(cavity=_find("cavity", cavity), eses=_find("eses", eses),
                 pdb2pqr=_find("pdb2pqr", pdb2pqr) if need_pdb2pqr else [])


def doctor() -> str:
    """Human-readable report of which tools were found."""
    lines = []
    for key, label in (("cavity", "CavityDetection"), ("eses", "ESES MS_Intersection"), ("pdb2pqr", "pdb2pqr")):
        try:
            cmd = _find(key)
            extra = ""
            if key == "pdb2pqr":
                try:
                    out = subprocess.run(cmd + ["--version"], capture_output=True, text=True, timeout=60)
                    extra = "  (" + (out.stdout or out.stderr).strip().splitlines()[0] + ")"
                except Exception:
                    pass
            lines.append(f"OK       {label:22s} {' '.join(cmd)}{extra}")
        except ToolNotFound as e:
            lines.append(f"MISSING  {label:22s} {e}")
    for mod in ("numpy", "scipy", "pandas"):
        ok = importlib.util.find_spec(mod) is not None
        lines.append(f"{'OK' if ok else 'optional' if mod == 'pandas' else 'MISSING':8s} python: {mod}")
    return "\n".join(lines)
