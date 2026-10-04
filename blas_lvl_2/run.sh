#!/bin/bash
# Build and run the BLAS-2 abstraction benchmark per compiler/flag set.
# -ffp-contract=off everywhere: FMA contraction would otherwise differ between
# a fused loop and an expression the compiler evaluates in separate passes,
# which is an arithmetic difference, not an abstraction penalty.
# Kernels and driver are compiled separately, no LTO.
set -e
mkdir -p build results
configs=(
  "gfortran|gfortran|-O2"
  "gfortran|gfortran|-O3 -march=native"
  "flang19|flang-new-19|-O2"
  "flang19|flang-new-19|-O3 -march=native"
)
for c in "${configs[@]}"; do
  IFS='|' read -r tag fc flags <<< "$c"
  flags="$flags -ffp-contract=off"
  id="${tag}_$(echo "$flags" | sed 's/-ffp-contract=off//' | tr -d ' =-')"
  d=build/$id; mkdir -p $d
  objs="$d/f77.o $d/modern.o"
  $fc $flags -c blas2_f77.f -o $d/f77.o
  $fc $flags -c blas2_modern.f90 -J $d -o $d/modern.o
  $fc $flags -S blas2_f77.f -o $d/f77.s
  $fc $flags -S blas2_modern.f90 -J $d -o $d/modern.s
  pdt=""
  if $fc $flags -c blas2_pdt.f90 -J $d -o $d/pdt.o 2>/dev/null; then
    pdt="-DHAVE_PDT"; objs="$objs $d/pdt.o"
    $fc $flags -S blas2_pdt.f90 -J $d -o $d/pdt.s
  fi
  $fc -O2 $pdt -c bench2.F90 -I $d -J $d -o $d/bench2.o
  $fc $d/bench2.o $objs -o $d/bench2
  echo "== $tag $flags ${pdt:-(no PDT support)}"
  { echo "# compiler: $tag"; echo "# flags: $flags"; [ -z "$pdt" ] && echo "# pdt: unsupported";
    taskset -c 1 $d/bench2; } > results/$id.txt
  grep '^#' results/$id.txt | grep -v -E 'compiler|flags' || true
done
