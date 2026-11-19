"""Result objects and parsers for the files the C++ programs write."""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# ------------------------------------------------------------------ parsers
def read_pqr(path: str, residue_order: Optional[Sequence[Tuple[str, str, str]]] = None) -> Dict[str, np.ndarray]:
    """Atoms of a PQR file in file order (= .xyzr order = CavityDetection atom indices).

    PDB2PQR 3.x drops chain IDs; pass residue_order [(chain, resnum, icode), ...] of the
    input PDB to restore them."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            head = line[:30].split()
            resnum = head[-1]
            name = line[12:16].strip()
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                f = line[54:].split()
                q, r = float(f[0]), float(f[1])
            except (ValueError, IndexError):
                f = line.split()
                xyz, q, r = (float(f[-5]), float(f[-4]), float(f[-3])), float(f[-2]), float(f[-1])
            chain = head[-2] if len(head) >= 6 and len(head[-2]) == 1 else ""   # pdb2pqr 3.x omits it
            icode = resnum[-1] if resnum[-1].isalpha() else ""
            rows.append((chain, re.sub(r"[A-Za-z]$", "", resnum), icode, line[17:20].strip(), name, xyz, q, r))
    if not rows:
        raise ValueError(f"no atoms in {path}")
    chain = [r[0] for r in rows]
    if residue_order:
        k, prev = -1, None
        for i, r in enumerate(rows):
            cur = (r[1], r[2])
            if cur != prev:
                k += 1
                prev = cur
            if 0 <= k < len(residue_order) and tuple(residue_order[k][1:]) == cur:
                chain[i] = residue_order[k][0]
    names = [r[4] for r in rows]
    return dict(
        chain=np.array(chain, dtype=object),
        resnum=np.array([int(r[1]) for r in rows]),
        icode=np.array([r[2] for r in rows], dtype=object),
        resname=np.array([r[3] for r in rows], dtype=object),
        name=np.array(names, dtype=object),
        xyz=np.array([r[5] for r in rows], dtype=float),
        charge=np.array([r[6] for r in rows], dtype=float),
        radius=np.array([r[7] for r in rows], dtype=float),
        heavy=np.array([not (n.startswith("H") or (n[:1].isdigit() and n[1:2] == "H")) for n in names]),
    )


def write_xyzr(xyz: np.ndarray, radius: np.ndarray, path: str) -> None:
    with open(path, "w") as out:
        for (x, y, z), r in zip(np.asarray(xyz, float), np.asarray(radius, float)):
            out.write(f"{x:.3f} {y:.3f} {z:.3f} {r:.4f}\n")


_POCKET_RE = re.compile(r"Pocket (\d+):\s*Area:\s*([\d.eE+-]+)\s*Volume:\s*([\d.eE+-]+)\s*"
                        r"Max Depth:\s*([\d.eE+-]+)\s*(?:Depth \(A\):\s*[\d.eE+-]+\s*)?Adjacent Atoms:\s*\n([\d ]*)")


def read_result_txt(path: str) -> List[dict]:
    """Parse CavityDetection's Result.txt -> [{index, area, volume, depth_steps, atoms}, ...]."""
    txt = open(path).read()
    return [dict(index=int(m.group(1)), area=float(m.group(2)), volume=float(m.group(3)),
                 depth_steps=float(m.group(4)), atoms=sorted({int(x) for x in m.group(5).split()}))
            for m in _POCKET_RE.finditer(txt)]


def voxel_from_log(text: str) -> Optional[float]:
    """Voxel edge (A) used by CavityDetection: cbrt(N / 1e5), N = first 'active voxels' count."""
    m = re.search(r"Total # active voxels of protein:\s*(\d+)", text)
    return (int(m.group(1)) / 100000.0) ** (1.0 / 3.0) if m else None


def read_obj(path: str) -> Tuple[np.ndarray, np.ndarray]:
    v, f = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                v.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                f.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.array(v, float).reshape(-1, 3), np.array(f, int).reshape(-1, 3)


# ------------------------------------------------------------------ objects
@dataclass
class Pocket:
    index: int                 # pocket number as in Result.txt
    atom_indices: np.ndarray   # 0-based indices into PocketResult.atoms (lining atoms)
    area: float                # A^2
    volume: float              # A^3
    depth_steps: float         # CavityDetection's "Max Depth" (deformation steps)
    depth: Optional[float]     # A  (= steps * 0.5 * voxel)
    residues: List[str]        # lining residues, e.g. "A:TYR26"
    centroid: np.ndarray       # mean of the lining atoms
    rank: int = 0              # 1 = largest volume

    @property
    def n_atoms(self) -> int:
        return len(self.atom_indices)

    def summary(self) -> dict:
        return dict(pocket=self.index, rank=self.rank, volume_A3=round(self.volume, 3), area_A2=round(self.area, 3),
                    depth_A=None if self.depth is None else round(self.depth, 3), depth_steps=self.depth_steps,
                    n_atoms=self.n_atoms, n_residues=len(self.residues),
                    centroid=[round(float(c), 3) for c in self.centroid], residues=self.residues)

    def __repr__(self) -> str:
        d = f"{self.depth:.2f} A" if self.depth is not None else f"{self.depth_steps:g} steps"
        return (f"Pocket({self.index}: volume {self.volume:.1f} A^3, area {self.area:.1f} A^2, depth {d}, "
                f"{len(self.residues)} residues)")


@dataclass
class PocketResult:
    pockets: List[Pocket]
    atoms: Dict[str, np.ndarray]           # receptor atoms (PQR order): chain, resnum, resname, name, xyz, ...
    voxel_size: Optional[float]
    settings: dict
    receptor_chains: List[str] = field(default_factory=list)
    peptide_chains: List[str] = field(default_factory=list)
    peptide_xyz: Optional[np.ndarray] = None   # peptide heavy atoms from the input file (if any)
    notes: List[str] = field(default_factory=list)
    log: str = ""
    workdir: Optional[str] = None              # kept only with keep_files=True
    mesh: Optional[Tuple[np.ndarray, np.ndarray]] = None   # ESES surface (vertices, faces) if requested

    # -------------------------------------------------- access
    def __len__(self):
        return len(self.pockets)

    def __iter__(self):
        return iter(self.pockets)

    def __getitem__(self, i) -> Pocket:
        return self.pockets[i]

    def __repr__(self) -> str:
        lines = [f"PocketResult: {len(self.pockets)} pocket(s), receptor chains {''.join(self.receptor_chains)}, "
                 f"{len(self.atoms['xyz'])} atoms, voxel {self.voxel_size or float('nan'):.3f} A"]
        lines += [f"  {p!r}" for p in self.by_volume()]
        return "\n".join(lines)

    def by_volume(self) -> List[Pocket]:
        return sorted(self.pockets, key=lambda p: -p.volume)

    def residue_label(self, i: int) -> str:
        a = self.atoms
        return f"{a['chain'][i]}:{a['resname'][i]}{a['resnum'][i]}{a['icode'][i]}"

    def pocket_atoms(self, pocket: "Pocket|int") -> Dict[str, np.ndarray]:
        p = self.pockets[pocket] if isinstance(pocket, int) else pocket
        return {k: v[p.atom_indices] for k, v in self.atoms.items()}

    # -------------------------------------------------- export
    def to_records(self) -> List[dict]:
        return [p.summary() for p in self.pockets]

    def to_dataframe(self):
        import pandas as pd
        df = pd.DataFrame(self.to_records())
        df["residues"] = df["residues"].apply(" ".join)
        return df.set_index("pocket")

    def atoms_dataframe(self):
        import pandas as pd
        a = self.atoms
        df = pd.DataFrame({k: list(v) for k, v in a.items() if k != "xyz"})
        df[["x", "y", "z"]] = a["xyz"]
        df["pocket"] = -1
        for p in self.pockets:
            df.loc[p.atom_indices, "pocket"] = p.index
        return df

    def to_json(self, path: Optional[str] = None) -> str:
        d = dict(settings=self.settings, voxel_size=self.voxel_size, receptor_chains=self.receptor_chains,
                 peptide_chains=self.peptide_chains, notes=self.notes,
                 pockets=[dict(p.summary(), atom_indices=[int(i) for i in p.atom_indices]) for p in self.pockets])
        s = json.dumps(d, indent=1)
        if path:
            with open(path, "w") as fh:
                fh.write(s)
        return s

    def to_csv(self, path: str) -> None:
        import csv
        recs = self.to_records()
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["pocket", "rank", "volume_A3", "area_A2", "depth_A", "depth_steps", "n_atoms",
                        "n_residues", "cx", "cy", "cz", "residues"])
            for r in recs:
                w.writerow([r["pocket"], r["rank"], r["volume_A3"], r["area_A2"], r["depth_A"], r["depth_steps"],
                            r["n_atoms"], r["n_residues"], *r["centroid"], " ".join(r["residues"])])

    def to_pdb(self, path: str) -> None:
        """Receptor PDB with the pocket number (+1, 0 = none) in the B-factor column."""
        a = self.atoms
        b = np.zeros(len(a["xyz"]))
        for p in self.pockets:
            b[p.atom_indices] = p.index + 1
        with open(path, "w") as out:
            for i in range(len(b)):
                x, y, z = a["xyz"][i]
                nm = a["name"][i]
                nm = f" {nm:<3}" if len(nm) < 4 else nm
                out.write(f"ATOM  {(i + 1) % 100000:5d} {nm} {a['resname'][i]:>3} {a['chain'][i] or ' '}"
                          f"{a['resnum'][i]:4d}{a['icode'][i] or ' '}   {x:8.3f}{y:8.3f}{z:8.3f}"
                          f"{1.0:6.2f}{b[i]:6.2f}\n")
            out.write("END\n")

    def to_pymol(self, path: str, receptor_pdb: str = "pockets.pdb") -> None:
        """PyMOL script colouring each pocket's lining atoms (load after the PDB from to_pdb)."""
        colors = ["red", "blue", "green", "orange", "magenta", "cyan", "yellow", "salmon", "purple", "lime"]
        name = os.path.splitext(os.path.basename(receptor_pdb))[0]
        lines = [f"load {receptor_pdb}", "hide everything", "show surface", "color grey80"]
        for k, p in enumerate(self.by_volume()):
            ids = "+".join(str(i + 1) for i in p.atom_indices)
            lines += [f"select pocket{p.index}, {name} and id {ids}", f"color {colors[k % len(colors)]}, pocket{p.index}"]
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    def save(self, outdir: str, prefix: str = "pockets") -> str:
        os.makedirs(outdir, exist_ok=True)
        self.to_json(os.path.join(outdir, f"{prefix}.json"))
        self.to_csv(os.path.join(outdir, f"{prefix}.csv"))
        self.to_pdb(os.path.join(outdir, f"{prefix}.pdb"))
        self.to_pymol(os.path.join(outdir, f"{prefix}.pml"), f"{prefix}.pdb")
        with open(os.path.join(outdir, f"{prefix}.log"), "w") as fh:
            fh.write(self.log)
        return outdir

    # -------------------------------------------------- evaluation against a known ligand
    def evaluate(self, ligand_xyz: Optional[np.ndarray] = None, cutoff: float = 4.5) -> List[dict]:
        """Compare every pocket with a bound ligand/peptide (heavy-atom coordinates, same frame).

        Defaults to the peptide chains found in the input file.  Binding residues = receptor
        residues with a heavy atom within `cutoff` A of the ligand.  Returns one dict per pocket
        (best residue F1 first): residue recall/precision/F1, atom recall/precision, DCC (A,
        pocket lining centroid to ligand centroid)."""
        lig = self.peptide_xyz if ligand_xyz is None else np.asarray(ligand_xyz, float).reshape(-1, 3)
        if lig is None or len(lig) == 0:
            raise ValueError("no ligand coordinates (input had no peptide chain; pass ligand_xyz)")
        from scipy.spatial import cKDTree
        a = self.atoms
        dist, _ = cKDTree(lig).query(a["xyz"])
        contact = set(np.where(a["heavy"] & (dist <= cutoff))[0].tolist())
        reskey = [self.residue_label(i) for i in range(len(dist))]
        binding = {reskey[i] for i in contact}
        lc = lig.mean(axis=0)
        rows = []
        for p in self.pockets:
            idx = p.atom_indices.tolist()
            lh = {i for i in idx if a["heavy"][i]}
            lres = {reskey[i] for i in idx}
            ov, rov = len(lh & contact), len(lres & binding)
            rr = rov / len(binding) if binding else 0.0
            rp = rov / len(lres) if lres else 0.0
            rows.append(dict(pocket=p.index, volume_A3=p.volume, residue_recall=rr, residue_precision=rp,
                             residue_f1=2 * rr * rp / (rr + rp) if rr + rp else 0.0,
                             atom_recall=ov / len(contact) if contact else 0.0,
                             atom_precision=ov / len(lh) if lh else 0.0,
                             dcc=float(np.linalg.norm(p.centroid - lc)),
                             shared_residues=sorted(lres & binding, key=_res_sort)))
        rows.sort(key=lambda r: (-r["residue_f1"], r["dcc"]))
        self.binding_residues = sorted(binding, key=_res_sort)
        return rows


def _res_sort(label: str):
    ch, rest = label.split(":", 1)
    m = re.search(r"(-?\d+)", rest)
    return (ch, int(m.group(1)) if m else 0, rest)
