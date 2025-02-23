#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -ex

CARGS=${CARGS:-""}

clargs(){
  # enable gpu if support is considered available
  set -e
  ARGS=""
  [ -d /opt/rocm/bin ] && ARGS="res/mkn/hip_mpi"
  [ -z "$ARGS" ] && which clang 2>&1 > /dev/null && clang -v 2>&1 | grep -q "Found CUDA" && ARGS="res/mkn/clang_cuda"
  [ -n "$ARGS" ] && ARGS="-P mkn.base=gpu_ -x $ARGS"
  echo "$ARGS"
}

# FILE="tests/core/data/tiles/test_tile.cpp"
# FILE="tests/core/data/particles/test_particles_selecting.cpp"
# FILE="tests/core/data/tiles/test_tile_set_mapper.cpp"
# FILE="tests/core/numerics/ampere/test_tileset_ampere.cpp"
# FILE="tests/core/numerics/faraday/test_tileset_faraday.cpp"
# FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
# FILE="tests/amr/data/particles/copy/test_particledata_copyNd.cpp"
# FILE="tests/amr/data/particles/stream_pack/test_main.cpp"
# FILE="tests/amr/data/particles/copy_overlap/test_particledata_copy_periodicNd.cpp"
# FILE="tests/amr/data/particles/test_particles_data.cpp"
TEST="$(clargs) -M ${FILE}"

set -x

(
  [ ! -f "bin/core/libphare_core.so" ] && (
      mkn clean build -qgp test_core -d google.test,+ -a "-fPIC" ${CARGS}
      mkn clean build -tqgp core -a "-fPIC" -x res/mkn/mpi ${CARGS}
  )
  # [ ! -f "bin/amr/libphare_amr.so" ] && (
  #     mkn clean build -qgp amr -a "-fPIC" -x res/mkn/mpi ${CARGS}
  # )

  mkn -p test_core ${TEST} ${CARGS} "$@"
  # mkn -p test_amr ${TEST} ${CARGS} "$@"

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
