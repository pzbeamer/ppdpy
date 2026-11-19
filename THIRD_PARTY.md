# Third-party software

ppdpy is a thin wrapper. The scientific methods and the programs that do the work are by
other authors; please cite them (see `CITATION.cff` and the README).

| Component | Authors | Source | Used as |
|---|---|---|---|
| **PPD** — `CavityDetection` (convex-hull surface evolution + Reeb graph pocket detection) | Rundong Zhao, Zixuan Cang, Yiying Tong, Guo-Wei Wei | https://github.com/rdzhao/PPD @ `4f90c45` (2018-05-29) | cloned and compiled by `ppd build-backend`; not redistributed |
| **ESES** — `MS_Intersection` (Eulerian solvent excluded surface) | Beibei Liu, Bao Wang, Rundong Zhao, Yiying Tong, Guo-Wei Wei | https://github.com/WeilabMSU/ESES @ `0b4cd51` (2020-08-04) | cloned and compiled by `ppd build-backend`; not redistributed |
| **PDB2PQR** — protonation and CHARMM radii/charges | Dolinsky, Nielsen, McCammon, Baker; Jurrus et al. | https://github.com/Electrostatics/pdb2pqr (BSD-3-Clause) | separate install (`pip install pdb2pqr`), called as a program |

Neither PPD nor ESES includes a license file upstream, so this repository contains none of
their source code. `ppdpy/backend/patches/` holds only small diffs written for this project,
applied to the fresh clones at build time:

- `ppd-0001-removepocket-return.patch` — adds a missing `return` at the end of
  `CavityDetection::RemovePocket`. Without it, recent GCC (13+) builds crash with
  "Illegal instruction" (falling off the end of a non-void function is undefined behaviour).
- `ppd-0002-portable-cmake.patch` — rewrites `CMakeLists.txt` for current toolchains on
  Linux and macOS (C++17, imported targets for TBB/Boost/CGAL/OpenVDB, optional OpenMP,
  `-fpermissive` for GCC only, `-ffp-contract=off` for clang so arm64 matches x86) and fixes the
  TBB version lookup in `cmake_modules/FindTBB.cmake` for oneTBB. The build script also removes
  the bundled `cmake_modules/FindBoost.cmake`, which fails with current Boost.
- `eses-0001-portable-cmake.patch` — finds `libomp` from conda or Homebrew on macOS instead of
  hard-coded `/usr/local` paths, Release build by default, `-ffp-contract=off` for clang.

None of the patches change the algorithms. Results were checked to be identical to builds
of the unpatched code on Linux (GCC) for the test structures.

The example structure `examples/1KNE.pdb` is from the Protein Data Bank (https://www.rcsb.org/structure/1KNE),
distributed by wwPDB under CC0 1.0.
