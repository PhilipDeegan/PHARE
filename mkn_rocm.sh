#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
FILE="tests/core/numerics/faraday/test_tileset_faraday.cpp"
FILE="tests/core/numerics/ampere/test_tileset_ampere.cpp"
FILE="tests/core/data/tiles/test_tile.cpp"
FILE="tests/amr/data/field/test_fields_schedules.cpp"
FILE="tests/amr/level_initializer/test_hybrid_level_initializer.cpp"
FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
CARGS=${CARGS:-""}
# CI (.github/workflows/build.yml) points these at generated profiles for its GPU (simulated)
MKN_HIP=${MKN_HIP:-"res/mkn/hip"} # res/mkn/hip_thrust: with rocThrust
MKN_HIP_MPI=${MKN_HIP_MPI:-"res/mkn/hip_mpi"}
MKN_HIP_RUN=${MKN_HIP_RUN:-"run"} # MKN_HIP_RUN=dbg for gdb

set -x
(
  [ ! -f "bin/core/libphare_core.a" ] && (
    set -e
    mkn clean build -Kqp test_core -d google.test,+ -Oa "-fPIC" ${CARGS}
    mkn clean build -Kqp core -a "-fPIC" -x res/mkn/mpi -Og 0
  )

  [ ! -f "bin/amr/libphare_amr.a" ] && (
    set -e
    mkn clean build -Kqdp test_amr -a "-fPIC" -x res/mkn/mpi -Og 0
  )

  # if/else, not `[[ ]] && ( ) || true`: a false test is no failure, a failed build/run still is
  if [[ $FILE == tests/core/* ]]; then
    mkn clean build -p test_core -M "${FILE}" ${CARGS} -w mkn.gpu -x ${MKN_HIP} ${MKN_HIP_RUN} "$@"
  else
    mkn clean build -qp test_diagnostics -M "${FILE}" ${CARGS} -w mkn.gpu -x ${MKN_HIP_MPI} run "$@"
    echo "ok"
  fi

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

# stop if not CI
[ -z "$CI" ] && exit 0 # continue for more

(
  set -xe
  mkn -tp core_tests ${CARGS} -w mkn.gpu -x ${MKN_HIP} clean build test run "$@"
  mkn -tp more_tests ${CARGS} -w mkn.gpu -x ${MKN_HIP_MPI} clean build test run "$@"
)
