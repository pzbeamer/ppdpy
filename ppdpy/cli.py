"""Command line:  ppd detect FILE.pdb [options]  |  ppd doctor  |  ppd build-backend"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from . import __version__
from .pipeline import PipelineError, detect_pockets
from .tools import ToolNotFound, backend_home, doctor


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="ppd", description="PPD pocket detection (pdb2pqr -> ESES -> CavityDetection)")
    ap.add_argument("--version", action="version", version=f"ppdpy {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("doctor", help="show which external programs were found")
    b = sub.add_parser("build-backend", help="fetch, patch and compile CavityDetection (PPD) and ESES")
    b.add_argument("--prefix", help=f"install dir, binaries in PREFIX/bin (default {backend_home()})")
    b.add_argument("--src", help="keep the upstream sources and build tree here")
    b.add_argument("-j", "--jobs", type=int)
    d = sub.add_parser("detect", help="detect pockets in a PDB or PQR file")
    d.add_argument("structure", help=".pdb (prepared automatically) or .pqr file")
    d.add_argument("-o", "--out", help="output folder (default: <file stem>_pockets)")
    d.add_argument("--chains", help="receptor chains, e.g. A or AB (default: chains > --peptide-max-len residues)")
    d.add_argument("--peptide", help="peptide chains to remove (default: remaining short chains)")
    d.add_argument("--model", type=int, default=1, help="model in NMR files (default 1)")
    d.add_argument("--residues", help="keep only these receptor residues, e.g. 20-80 or 5-40,45-80")
    d.add_argument("--peptide-max-len", type=int, default=30)
    d.add_argument("--probe", type=float, default=1.4)
    d.add_argument("--resolution", type=float, default=0.4)
    d.add_argument("--margin", type=float, default=1.0)
    d.add_argument("--min-depth", type=float, default=4)
    d.add_argument("--min-area", type=float, default=80)
    d.add_argument("--min-volume", type=float, default=150)
    d.add_argument("--cutoff", type=float, default=4.5, help="contact cut-off for scoring against the peptide")
    d.add_argument("--keep-files", action="store_true", help="keep pqr/obj/xyzr/Result.txt in OUT/work")
    d.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args(argv)

    if a.cmd == "doctor":
        print(doctor())
        return 0
    if a.cmd == "build-backend":
        script = Path(__file__).with_name("backend") / "build_backend.sh"
        cmd = ["bash", str(script)]
        if a.prefix:
            cmd += ["--prefix", a.prefix]
        if a.src:
            cmd += ["--src", a.src]
        if a.jobs:
            cmd += ["--jobs", str(a.jobs)]
        return subprocess.call(cmd)

    out = a.out or os.path.splitext(os.path.basename(a.structure))[0] + "_pockets"
    try:
        res = detect_pockets(a.structure, chains=a.chains, peptide=a.peptide, model=a.model, residues=a.residues,
                             peptide_max_len=a.peptide_max_len, probe=a.probe, resolution=a.resolution,
                             margin=a.margin, min_depth=a.min_depth, min_area=a.min_area, min_volume=a.min_volume,
                             workdir=os.path.join(out, "work") if a.keep_files else None, verbose=not a.quiet)
    except (ToolNotFound, PipelineError, ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    res.save(out)
    print(res)
    for n in res.notes:
        print("note:", n)
    if res.peptide_xyz is not None:
        rows = res.evaluate(cutoff=a.cutoff)
        print(f"\nagainst peptide chain(s) {''.join(res.peptide_chains)} "
              f"({len(res.binding_residues)} binding residues within {a.cutoff} A):")
        for r in rows:
            print(f"  pocket {r['pocket']}: residue recall {r['residue_recall']:.2f}  precision "
                  f"{r['residue_precision']:.2f}  F1 {r['residue_f1']:.2f}  DCC {r['dcc']:.1f} A")
    print(f"\nwrote {out}/pockets.{{json,csv,pdb,pml,log}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
