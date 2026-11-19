"""PDB preparation: pick a model and chains, drop waters/ligands/hydrogens/altlocs,
convert modified residues of the receptor to their parent amino acid, and split off
short (peptide) chains.  Same rules as scripts/prep_complex.py."""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
}
BACKBONE = ["N", "CA", "C", "O", "OXT"]
SIDE = {
    "ALA": ["CB"], "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"], "ASN": ["CB", "CG", "OD1", "ND2"],
    "ASP": ["CB", "CG", "OD1", "OD2"], "CYS": ["CB", "SG"], "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"], "GLY": [], "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"], "LEU": ["CB", "CG", "CD1", "CD2"], "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "MET": ["CB", "CG", "SD", "CE"], "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CB", "CG", "CD"], "SER": ["CB", "OG"], "THR": ["CB", "OG1", "CG2"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], "VAL": ["CB", "CG1", "CG2"],
}
MODIFIED = {
    "MSE": "MET", "M3L": "LYS", "MLY": "LYS", "MLZ": "LYS", "ALY": "LYS", "KCX": "LYS", "LLP": "LYS",
    "SEP": "SER", "TPO": "THR", "PTR": "TYR", "CSO": "CYS", "CSD": "CYS", "OCS": "CYS", "CME": "CYS",
    "CSS": "CYS", "HYP": "PRO", "NLE": "LEU", "FME": "MET", "2MR": "ARG", "AGM": "ARG", "DA2": "ARG",
    "MEN": "ASN", "SAC": "SER", "CGU": "GLU", "PCA": "GLU",
}
RENAME = {("MSE", "SE"): ("SD", "S")}
WATERS = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}

ResidueSpec = Union[None, str, Sequence[Tuple[int, int]], Sequence[int], range]


def parse_ranges(spec: ResidueSpec) -> List[Tuple[int, int]]:
    """'20-80' / '5-40,45-80' / range(20, 81) / [(20, 80)] / [20, 21, 22] -> [(a, b), ...]"""
    if spec is None or spec == "":
        return []
    if isinstance(spec, str):
        out = []
        for part in filter(None, spec.replace(" ", "").split(",")):
            a, _, b = part.partition("-")
            out.append((int(a), int(b or a)))
        return out
    if isinstance(spec, range):
        return [(spec.start, spec.stop - 1)]
    out = []
    for item in spec:
        if isinstance(item, (tuple, list)):
            out.append((int(item[0]), int(item[1])))
        else:
            out.append((int(item), int(item)))
    return out


def _read_lines(path_or_text: str, model: int):
    """ATOM/HETATM lines of one model, first altloc, no hydrogens."""
    if "\n" in path_or_text:
        lines = path_or_text.splitlines()
    else:
        with open(path_or_text) as fh:
            lines = fh.read().splitlines()
    current, n_models, out = None, 0, []
    for line in lines:
        rec = line[:6]
        if rec == "MODEL ":
            n_models += 1
            current = n_models
            continue
        if rec == "ENDMDL":
            current = None
            continue
        if rec not in ("ATOM  ", "HETATM"):
            continue
        if current is not None and current != model:
            continue
        if current is None and n_models and model != 1:
            continue
        line = line.ljust(80)
        if line[16] not in (" ", "A"):
            continue
        element = line[76:78].strip()
        name = line[12:16].strip()
        if element in ("H", "D") or (not element and (name[:1] in ("H", "D") or name[:2] in ("1H", "2H", "3H"))):
            continue
        out.append(line[:16] + " " + line[17:])
    return out, max(n_models, 1)


@dataclass
class PreparedStructure:
    receptor_pdb: str                      # PDB text of the receptor (input to pdb2pqr)
    peptide_pdb: str                       # PDB text of the peptide chains ("" if none)
    receptor_chains: List[str]
    peptide_chains: List[str]
    peptide_xyz: np.ndarray                # (M, 3) heavy atoms of the peptide chains
    residue_order: List[Tuple[str, str, str]]   # (chain, resnum, icode) of receptor residues, file order
    notes: List[str] = field(default_factory=list)
    n_models: int = 1


def prepare_structure(pdb: str, chains: Optional[Iterable[str]] = None, peptide: Optional[Iterable[str]] = None,
                      model: int = 1, residues: ResidueSpec = None, peptide_max_len: int = 30) -> PreparedStructure:
    """Split a PDB file (path or text) into receptor and peptide.

    chains   receptor chain IDs (default: every protein chain longer than peptide_max_len,
             or all chains if none is that long)
    peptide  peptide chain IDs (default: the other short protein chains)
    residues keep only these receptor residue numbers ('20-80', [(5, 40), (45, 80)], range(...))
    """
    lines, n_models = _read_lines(pdb, model)
    if model > n_models:
        raise ValueError(f"model {model} requested but the file has {n_models} model(s)")
    res_lines: "OrderedDict[tuple, list]" = OrderedDict()
    for line in lines:
        res_lines.setdefault((line[21], line[22:26].strip(), line[26]), []).append(line)

    chain_res: "OrderedDict[str, OrderedDict]" = OrderedDict()
    dropped: Dict[str, int] = {}
    for key, ls in res_lines.items():
        resn = ls[0][17:20].strip()
        names = {l[12:16].strip() for l in ls}
        if resn in THREE_TO_ONE and ls[0].startswith("ATOM"):
            kind = "std"
        elif {"N", "CA", "C"} <= names and (resn in MODIFIED or resn in THREE_TO_ONE):
            kind = "mod"
        elif {"N", "CA", "C"} <= names and key[0] in chain_res:
            kind = "mod?"
        else:
            dropped[resn] = dropped.get(resn, 0) + len(ls)
            continue
        chain_res.setdefault(key[0], OrderedDict())[key] = (resn, ls, kind)
    if not chain_res:
        raise ValueError("no protein residues found")

    lengths = {c: len(r) for c, r in chain_res.items()}
    rec = list(chains) if chains else []
    pep = list(peptide) if peptide else []
    if not rec:
        rec = [c for c in chain_res if c not in pep and lengths[c] > peptide_max_len] or \
              [c for c in chain_res if c not in pep]
    if peptide is None:
        pep = [c for c in chain_res if c not in rec and lengths[c] <= peptide_max_len]
    missing = [c for c in rec + pep if c not in chain_res]
    if missing:
        raise ValueError(f"chain(s) {''.join(missing)} not found; protein chains: {''.join(chain_res)}")

    ranges = parse_ranges(residues)
    notes: List[str] = []
    rec_lines: List[str] = []
    order: List[Tuple[str, str, str]] = []
    for ch in rec:
        for (c, num, ic), (resn, ls, kind) in chain_res[ch].items():
            if ranges and not any(a <= int(num) <= b for a, b in ranges):
                continue
            if kind == "std":
                rec_lines += ls
                order.append((c, num, ic.strip()))
                continue
            parent = MODIFIED.get(resn) or (resn if resn in THREE_TO_ONE else None)
            if parent is None:
                notes.append(f"receptor {resn}{num}{ic.strip()}:{ch} dropped (unknown modified residue)")
                continue
            keep = set(BACKBONE + SIDE.get(parent, []))
            for l in ls:
                nm = l[12:16].strip()
                new, el = RENAME.get((resn, nm), (nm, None))
                if new not in keep:
                    continue
                l = "ATOM  " + l[6:12] + f" {new:<3}" + l[16:17] + f"{parent:>3}" + l[20:]
                if el:
                    l = l[:76] + f"{el:>2}" + l[78:]
                rec_lines.append(l)
            order.append((c, num, ic.strip()))
            notes.append(f"receptor {resn}{num}{ic.strip()}:{ch} converted to {parent}")
        rec_lines.append("TER")
    if not order:
        raise ValueError("no receptor residues left after chain/residue selection")

    pep_lines: List[str] = []
    for ch in pep:
        for _, (resn, ls, kind) in chain_res[ch].items():
            pep_lines += ls
        pep_lines.append("TER")
    pxyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in pep_lines if l != "TER"],
                    dtype=float).reshape(-1, 3)

    other = [c for c in chain_res if c not in rec and c not in pep]
    if other:
        notes.append(f"ignored chains: {''.join(other)}")
    if dropped:
        notes.append("dropped: " + ", ".join(f"{k} ({v} atoms)" for k, v in dropped.items()))
    if n_models > 1:
        notes.append(f"model {model} of {n_models}")

    def text(ls):
        return "".join(l.rstrip() + "\n" for l in ls) + "END\n" if ls else ""

    return PreparedStructure(receptor_pdb=text(rec_lines), peptide_pdb=text(pep_lines),
                             receptor_chains=rec, peptide_chains=pep, peptide_xyz=pxyz,
                             residue_order=order, notes=notes, n_models=n_models)


def read_coordinates(pdb: str, chains: Optional[Iterable[str]] = None, model: int = 1) -> np.ndarray:
    """(N, 3) heavy-atom coordinates (waters excluded) of a PDB file (path or text),
    e.g. a ligand/peptide for evaluate()."""
    lines, _ = _read_lines(pdb, model)
    chains = set(chains) if chains else None
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines
           if (chains is None or l[21] in chains) and l[17:20].strip() not in WATERS]
    return np.array(xyz, dtype=float).reshape(-1, 3)
