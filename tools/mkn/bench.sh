#!/usr/bin/env bash

set -e

FILES=(
# "tools/bench/core/data/particles/sorting/bench_gpu.cpp"
  "tools/bench/core/numerics/interpolator/bench_interpolator_gpu.cpp"
)

CARG="-DPHARE_FORCE_DEBUG_DO=1"

# default cuda
XFIL="-x res/mkn/clang_cuda.yaml"

# enable hipcc if found
( which hipcc || [ -f /opt/rocm/bin/hipcc ] ) && XFIL="-x res/mkn/hip.yaml"

WITH="-w mkn.gpu"
export MKN_DBG="ncu -o profile"


for FILE in ${FILES[@]}; do
    echo $FILE
    mkn clean build -M "${FILE}" -a "${CARG}" -p bench ${DOPT} ${WITH} ${XFIL} run -Og 0 $@
    #-- --benchmark_repetitions=25 --benchmark_display_aggregates_only=true
done
