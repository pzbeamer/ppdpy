# ppdpy

**A Python wrapper for [PPD](https://github.com/rdzhao/PPD), protein pocket detection by
convex hull surface evolution and Reeb graphs.** You give it a PDB file and get the pockets back
as Python objects.

> **Credit.** The pocket-detection method and the `CavityDetection` program are the work of
> **Rundong Zhao, Zixuan Cang, Yiying Tong and Guo-Wei Wei** (Michigan State University):
>
> R. Zhao, Z. Cang, Y. Tong, G.-W. Wei. *Protein pocket detection via convex hull surface evolution
> and associated Reeb graph.* **Bioinformatics** 34(17):i830–i837 (2018).
> [doi:10.1093/bioinformatics/bty598](https://doi.org/10.1093/bioinformatics/bty598) ·
> code: [github.com/rdzhao/PPD](https://github.com/rdzhao/PPD)
>
> The molecular surface comes from **ESES** by B. Liu, B. Wang, R. Zhao, Y. Tong and G.-W. Wei
> ([J. Comput. Chem. 38:446–466, 2017](https://doi.org/10.1002/jcc.24682);
> [github.com/WeilabMSU/ESES](https://github.com/WeilabMSU/ESES)). Charges and radii come from
> **[PDB2PQR](https://github.com/Electrostatics/pdb2pqr)**.
>
> ppdpy only wraps these programs. **If you use it, please cite the PPD paper**
> (see `CITATION.cff`). This project is not affiliated with or endorsed by the original authors.

```
PDB ─prepare─► receptor.pdb ─PDB2PQR (CHARMM)─► .pqr ─ESES─► surface .obj ─PPD CavityDetection─► pockets
```

ppdpy does not re-implement anything. It runs the original C++ programs, so the pockets are the
same as running PPD by hand. It adds:

- structure preparation: NMR model, chains, residue ranges, removal of waters, ligands and
  alternate locations, and modified residues
- automatic file handling in a temporary directory
- parsing of results, including depth in Å and pdb2pqr's dropped chain IDs
- scoring against a bound ligand or peptide
- exports to JSON, CSV, PDB and PyMOL

## Installation

PPD and ESES are compiled from their own repositories when you build the backend. They are
**not** included here (see [THIRD_PARTY.md](THIRD_PARTY.md)).

### 1. Dependencies

The easiest way to get the dependencies is conda, on Linux or macOS:

```bash
git clone https://github.com/<you>/ppdpy.git && cd ppdpy
conda env create -f environment.yml
conda activate ppdpy
```

Alternatively, use system packages:

- **Ubuntu/Debian:** `sudo apt install git cmake clang libomp-dev libopenvdb-dev libtbb-dev libboost-dev libcgal-dev libgmp-dev libmpfr-dev`
- **macOS (Homebrew):** `brew install cmake libomp openvdb tbb boost cgal`. The conda route is better tested on macOS.

### 2. Install ppdpy and build the backend

```bash
pip install -e ".[all]"      # ppdpy + pdb2pqr + pandas
ppd build-backend            # clones PPD and ESES at pinned commits, patches, compiles (a few minutes)
ppd doctor                   # checks that CavityDetection, MS_Intersection and pdb2pqr are found
```

The binaries are installed to `~/.local/share/ppdpy/bin`. To change that location, set `PPDPY_HOME`.
If you already have compiled binaries, you can skip the build and point ppdpy at them with
`PPD_CAVITYDETECTION`, `PPD_ESES` and `PPD_PDB2PQR`, or put them on your `PATH`.

## Quick start

```python
import ppdpy

res = ppdpy.detect_pockets("examples/1KNE.pdb")   # receptor chain A; short chain P is treated as the peptide
print(res)
# PocketResult: 2 pocket(s), receptor chains A, 880 atoms, voxel 0.568 A
#   Pocket(1: volume 424.4 A^3, area 313.6 A^2, depth 2.84 A, 17 residues)
#   Pocket(0: volume 211.1 A^3, area 163.2 A^2, depth 1.99 A, 9 residues)

p = res.by_volume()[0]
p.volume, p.area, p.depth          # Å³, Å², Å
p.residues                         # ['A:GLU1', 'A:TYR2', 'A:ALA3', ...]
p.atom_indices, p.centroid         # indices into res.atoms, (3,) array
res.pocket_atoms(p)["xyz"]         # lining-atom coordinates

res.evaluate()                     # vs the peptide in the file: residue recall/precision/F1, DCC
res.save("1KNE_pockets")           # pockets.json / .csv / .pdb (pocket id in B-factor) / .pml / .log
res.to_dataframe()                 # pandas
```

An unbound structure, a single domain of an NMR ensemble, or custom settings:

```python
apo = ppdpy.detect_pockets("1AP0.pdb", model=1, residues="20-80",
                           probe=1.4, resolution=0.4, min_depth=4, min_area=80, min_volume=150)
apo.evaluate(ligand_xyz=peptide_xyz)   # any (N, 3) array in the same coordinate frame
```

### `detect_pockets` options

| argument | default | meaning |
|---|---|---|
| `chains`, `peptide` | auto | Receptor chains, and chains removed as the bound peptide. By default, chains of `peptide_max_len` (30) residues or fewer count as peptide. |
| `model`, `residues` | 1, all | NMR model to use, and receptor residues to keep (`"20-80"`, `"5-40,45-80"`) |
| `probe`, `resolution`, `margin` | 1.4, 0.4, 1.0 Å | ESES surface settings |
| `min_depth`, `min_area`, `min_volume` | 4, 80 Å², 150 Å³ | PPD thresholds. `min_depth` is in evolution steps and changes the segmentation itself. |
| `workdir`, `keep_files` | temp dir | Keep the `.pqr`, `.xyzr`, `marching_cubes.obj` and `Result.txt` files |
| `load_mesh` | False | Set `res.mesh` to the ESES surface as `(vertices, faces)` |
| `timeout`, `verbose` | – | Time limit in seconds for each program, and progress output |

You can also pass a `.pqr` file directly, which skips preparation and PDB2PQR.
The lower-level steps are available as `prepare_structure`, `run_pdb2pqr`, `run_eses` and
`run_cavity_detection`. `run_cavity_detection` accepts mesh and atom arrays as well as files.
Each run uses its own temporary directory, so you can run several in parallel.

## Command line

```bash
ppd detect examples/1KNE.pdb -o 1KNE_pockets
ppd detect 1AP0.pdb --residues 20-80 --probe 1.4 --keep-files
ppd doctor
ppd build-backend [--prefix DIR] [--jobs N]
```

## Notes

- **Depth.** Depth in Å is the `Max Depth` steps × 0.5 × voxel size. The voxel size is
  (active voxels / 10⁵)^(1/3), as PPD computes it, and is available as `res.voxel_size`.
- **Atom indices.** Indices follow the PQR atom order, which includes the hydrogens added by
  PDB2PQR. Chain IDs, which PDB2PQR 3.x drops, are restored from the input file.
- **Structure preparation.** ppdpy uses the first alternate location of each atom, drops
  hydrogens, waters and ligands, and converts modified receptor residues to their parent amino
  acid (MSE→MET, M3L→LYS, SEP→SER, …).
- **Exit codes.** `CavityDetection` exits with status 1 on success, and ppdpy accepts this.
- **Patches.** The build applies a few small patches to the upstream code (see THIRD_PARTY.md).
  These are build fixes for current compilers plus one bug fix: a missing `return` that makes
  GCC 13+ builds crash. They do not change the algorithm.

## Tests

```bash
pytest -q                 # fast tests; the end-to-end 1KNE test runs if the backend is built
```

## License

The wrapper code is under the MIT license (see [LICENSE](LICENSE)). PPD and ESES remain under their
authors' terms and are not distributed with this repository.
