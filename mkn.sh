#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

CARGS=${CARGS:-""}

clargs(){
  # enable gpu if support is considered available
  set -e
  ARGS=""
  [ -d /opt/rocm/bin ] && ARGS="res/mkn/hip_mpi"
  [ -z "$ARGS" ] && which clang 2>&1 > /dev/null && clang -v 2>&1 | grep -q "Found CUDA" && ARGS="res/mkn/clang_cuda"
  [ -n "$ARGS" ] && ARGS="-P mkn.base=gpu_ -x $ARGS"
  [ -z "$ARGS" ] && ARGS="-x res/mkn/mpi" # default
  echo "$ARGS"
}

# FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
# FILE="tests/amr/data/particles/copy/test_particledata_copyNd.cpp"
# FILE="tests/amr/data/particles/stream_pack/test_main.cpp"
# FILE="tests/amr/data/particles/copy_overlap/test_particledata_copy_periodicNd.cpp"
# FILE="tests/amr/data/particles/test_particles_data.cpp"
# FILE="tests/amr/data/particles/refine/test_particles_data_split.cpp"
# FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
FILE="tests/core/numerics/ion_updater/test_updater_pp_main.cpp"
TEST=" -M ${FILE} "
CARGS="${CARGS} $(clargs)"

(
  set -x
  [ ! -f "bin/core/libphare_core.a" ] && (
      mkn clean build -KOgqtp test_core -d google.test,+ -a "-fPIC" ${CARGS}
      mkn clean build -KOgqtp core -a "-fPIC" ${CARGS}
  )
  # [ ! -f "bin/amr/libphare_amr.a" ] && (
  #     mkn clean build -Kqtp amr -a "-fPIC" ${CARGS}
  # )

  mkn -p test_core ${TEST} ${CARGS} "$@"
  # mkn -p test_amr ${TEST} ${CARGS} "$@"

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

exit 0 # continue for more

(
  set -x
  arr=(
    tests/core/data/particles/test_particles_construction.cpp
    tests/core/data/particles/sorting/test_particle_sorting.cpp
    tests/core/data/tiles/test_tile.cpp
    tests/core/data/particles/test_particles_selecting.cpp
    tests/core/numerics/ampere/test_tileset_ampere.cpp
    tests/core/numerics/faraday/test_tileset_faraday.cpp
    tests/core/numerics/ohm/test_tileset_ohm.cpp
    tests/core/data/particles/test_particles_serialization.cpp
  )

  for i in ${!arr[@]}; do
      mkn -p test_core -M ${arr[i]} ${CARGS} clean build run -gO
  done
)
