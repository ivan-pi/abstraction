#!/bin/bash
# Build and run the benchmark for each compiler/flag configuration.
# Kernels and driver are compiled separately (no LTO) so nothing inlines.
set -e
mkdir -p build results
configs=(
  "gfortran|gfortran|-O2"
  "gfortran|gfortran|-O3 -march=native"
  "gfortran|gfortran|-O3 -march=native -ffast-math"
  "flang19|flang-new-19|-O2"
  "flang19|flang-new-19|-O3 -march=native"
  "flang19|flang-new-19|-O3 -march=native -ffast-math"
)
for c in "${configs[@]}"; do
  IFS='|' read -r tag fc flags <<< "$c"
  id="${tag}_$(echo "$flags" | tr -d ' =-' )"
  d=build/$id; mkdir -p $d
  $fc $flags -c blas1_loop.f   -o $d/loop.o
  $fc $flags -c blas1_array.f90 -o $d/array.o
  $fc $flags -S blas1_loop.f   -o $d/loop.s
  $fc $flags -S blas1_array.f90 -o $d/array.s
  $fc -O2 -c bench.f90 -J $d -o $d/bench.o
  $fc $d/bench.o $d/loop.o $d/array.o -o $d/bench
  echo "== $tag $flags"
  { echo "# compiler: $tag"; echo "# flags: $flags"; taskset -c 1 $d/bench; } > results/$id.txt
  grep '^#' results/$id.txt | grep -v -E 'compiler|flags'
done
