#!/usr/bin/env bash
# Build and run the CSR SpMV abstraction benchmark for each compiler at -O2
# and -O3 with native tuning, then write results/summary.txt.
#
# Works on Linux and macOS.
#
#   ./run.sh                              # gfortran and flang, if found
#   FCS="gfortran-15 flang" ./run.sh      # choose compilers explicitly
#   PIN="" ./run.sh                       # disable core pinning on Linux
#
# -ffp-contract=off everywhere so FMA fusion can't differ between variants.
# Kernels and driver compiled separately, no LTO.
set -euo pipefail
cd "$(dirname "$0")"

# Core pinning: taskset on Linux.  macOS has no user-level CPU pinning; a
# busy single-threaded process normally runs on a performance core anyway.
if [[ -z "${PIN+set}" ]]; then
  if command -v taskset >/dev/null 2>&1; then PIN="taskset -c 1"; else PIN=""; fi
fi

# Compilers: first gfortran and first flang found on PATH, unless FCS is set.
first_found() {
  local c
  for c in "$@"; do
    if command -v "$c" >/dev/null 2>&1; then echo "$c"; return; fi
  done
}
if [[ -z "${FCS:-}" ]]; then
  FCS="$(first_found gfortran gfortran-15 gfortran-14 gfortran-13) \
       $(first_found flang flang-new flang-new-21 flang-new-20 flang-new-19)"
fi
FCS=$(echo $FCS)   # squeeze whitespace
[[ -n "$FCS" ]] || { echo "no Fortran compiler found; set FCS" >&2; exit 1; }

# Native-tuning flag the compiler accepts: -march=native (x86, LLVM on
# arm64) or -mcpu=native (GCC on arm64); empty if neither works.
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf 'end\n' > "$tmp/t.f90"
native_flag() {
  local f
  for f in -march=native -mcpu=native; do
    if "$1" $f -c "$tmp/t.f90" -o "$tmp/t.o" >/dev/null 2>&1; then echo "$f"; return; fi
  done
}

mkdir -p build results
rm -f results/*.txt
for fc in $FCS; do
  native=$(native_flag "$fc")
  for opt in "-O2" "-O3 $native"; do
    flags="$opt -ffp-contract=off"
    id="$(basename "$fc")_$(echo "$opt" | tr -d ' =-')"
    d="build/$id"
    # Build from copies of the sources inside a fresh directory, so the
    # compiler can only find the module files it just wrote (a stray .mod
    # from another compiler next to the sources would otherwise be picked
    # up first, since the current directory is searched before -I paths).
    rm -rf "$d"
    mkdir -p "$d"
    cp spmv_f77.f spmv_modern.f90 bench_spmv.f90 "$d/"
    fcpath=$(command -v "$fc")
    ( cd "$d"
      "$fcpath" $flags -c spmv_f77.f -o f77.o
      "$fcpath" $flags -c spmv_modern.f90 -o modern.o
      "$fcpath" $flags -S spmv_f77.f -o f77.s
      "$fcpath" $flags -S spmv_modern.f90 -o modern.s
      "$fcpath" -O2 -c bench_spmv.f90 -o bench.o
      "$fcpath" bench.o f77.o modern.o -o bench_spmv )
    echo "== $fc $flags${PIN:+  (pinned: $PIN)}"
    { echo "# compiler: $(basename "$fc")"
      echo "# flags: $flags"
      $PIN "$d/bench_spmv"; } > "results/$id.txt"
    grep -E '^# (note|FAIL)' "results/$id.txt" || true
  done
done

python3 report.py > results/summary.txt
echo "wrote results/summary.txt"
