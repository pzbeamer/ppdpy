#!/usr/bin/env bash
# Fetch, patch and build the C++ backend of ppdpy:
#   CavityDetection  from PPD  (Zhao et al. 2018)  https://github.com/rdzhao/PPD
#   MS_Intersection  from ESES (Liu et al. 2017)   https://github.com/WeilabMSU/ESES
# The upstream sources are cloned at pinned commits; they are not part of this repository.
#
# Usage:  bash build_backend.sh [--prefix DIR] [--src DIR] [--jobs N] [--keep-src]
#   --prefix  install directory; binaries go to DIR/bin (default: ${PPDPY_HOME:-~/.local/share/ppdpy})
#   --src     where to clone/build (default: a temporary folder, removed afterwards)
# Environment: CC / CXX choose the compiler; CMAKE_ARGS adds CMake options.
# Needs: git, cmake >= 3.16, a C++17 compiler, OpenMP, OpenVDB, TBB, Boost (headers), CGAL.
set -euo pipefail

PPD_URL="https://github.com/rdzhao/PPD.git"
PPD_COMMIT="4f90c45c985bb21bd8d0ff7ce266a5c9dbd5536c"     # 2018-05-29
ESES_URL="https://github.com/WeilabMSU/ESES.git"
ESES_COMMIT="0b4cd5129702cb2f29a18e970e07e850bb24fc92"    # 2020-08-04

HERE="$(cd "$(dirname "$0")" && pwd)"
PATCHES="${HERE}/patches"
PREFIX="${PPDPY_HOME:-${HOME}/.local/share/ppdpy}"
SRC=""
KEEP=0
JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 2)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2 ;;
    --src) SRC="$2"; KEEP=1; shift 2 ;;
    --jobs|-j) JOBS="$2"; shift 2 ;;
    --keep-src) KEEP=1; shift ;;
    -h|--help) sed -n 2,12p "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

for t in git cmake; do
  command -v "$t" >/dev/null || { echo "error: $t not found" >&2; exit 1; }
done
if [[ -z "${SRC}" ]]; then
  SRC="$(mktemp -d "${TMPDIR:-/tmp}/ppdpy_build.XXXXXX")"
fi
mkdir -p "${SRC}" "${PREFIX}/bin"
cleanup() { [[ "${KEEP}" == "1" ]] || rm -rf "${SRC}"; }
trap cleanup EXIT

fetch() {   # url commit dir
  if [[ ! -d "$3/.git" ]]; then
    git clone --quiet "$1" "$3"
  fi
  git -C "$3" fetch --quiet origin "$2" 2>/dev/null || true
  git -C "$3" checkout --quiet --force "$2"
  git -C "$3" clean -qfdx
}

echo "==== fetching sources into ${SRC}"
fetch "${PPD_URL}" "${PPD_COMMIT}" "${SRC}/PPD"
fetch "${ESES_URL}" "${ESES_COMMIT}" "${SRC}/ESES"

echo "==== patching"
for p in "${PATCHES}"/ppd-*.patch; do
  echo "  PPD  $(basename "$p")"; git -C "${SRC}/PPD" apply --whitespace=nowarn "$p"
done
# upstream ships an old FindBoost.cmake that breaks with current Boost; CMake's own is used instead
rm -f "${SRC}/PPD/cmake_modules/FindBoost.cmake"
for p in "${PATCHES}"/eses-*.patch; do
  echo "  ESES $(basename "$p")"; git -C "${SRC}/ESES" apply --whitespace=nowarn "$p"
done

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  CMAKE_ARGS="-DCMAKE_PREFIX_PATH=${CONDA_PREFIX} ${CMAKE_ARGS:-}"
fi

echo "==== building ESES (MS_Intersection)"
cmake -S "${SRC}/ESES" -B "${SRC}/ESES/build" -DCMAKE_BUILD_TYPE=Release ${CMAKE_ARGS:-}
cmake --build "${SRC}/ESES/build" --parallel "${JOBS}"

echo "==== building PPD (CavityDetection)"
cmake -S "${SRC}/PPD" -B "${SRC}/PPD/build" -DCMAKE_BUILD_TYPE=Release ${CMAKE_ARGS:-}
cmake --build "${SRC}/PPD/build" --parallel "${JOBS}"

install -m 755 "${SRC}/ESES/build/MS_Intersection" "${PREFIX}/bin/MS_Intersection"
install -m 755 "${SRC}/PPD/build/CavityDetection" "${PREFIX}/bin/CavityDetection"
echo "==== installed"
echo "  ${PREFIX}/bin/CavityDetection"
echo "  ${PREFIX}/bin/MS_Intersection"
if [[ "${PREFIX}" != "${PPDPY_HOME:-${HOME}/.local/share/ppdpy}" ]]; then
  echo "note: non-default prefix; export PPDPY_HOME=${PREFIX} so ppdpy finds the binaries"
fi
