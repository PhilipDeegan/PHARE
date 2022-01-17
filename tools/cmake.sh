#!/usr/bin/env bash

# usage:
#  copy this file to tools/cmake.sh - and edit as you wish
#    tools/cmake.sh is ignored by git
#
set -ex
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root
exec 19>$CWD/.cmake.sh.cmd # set -x redirect
export BASH_XTRACEFD=19  # set -x redirect
THREADS=${THREADS:="11"}
BUILD_DIR=${BUILD_DIR:="$CWD/build"}
SAMRAI=${SAMRAI:="/mkn/r/llnl/samrai/master"} # "" = as subproject
FFF=("${BUILD_DIR}")
CMAKE_CONFIG="-DdevMode=ON -Dasan=OFF -Dbench=OFF -DwithCaliper=OFF -DtestMPI=OFF"
CMAKE_CXX_FLAGS="-g3 -DPHARE_DIAG_DOUBLES=1 -O0" #3 -march=native -mtune=native"
set -xe
time (
  date
  [ -n "$CLEAN" ] && (( $CLEAN == 1 )) && for f in ${FFF[@]}; do rm -rf $f; done
  [ ! -f "$CWD/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
  [ ! -d "$CWD/subprojects/cppdict/include" ] && git submodule update --init
  mkdir -p ${BUILD_DIR}
  [[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"
  ( export CC=gcc CXX=g++
    cd ${BUILD_DIR} && cmake $CWD ${SAMRAI} -G Ninja ${CMAKE_CONFIG} \
    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" ) # -DCMAKE_BUILD_TYPE=Release
  # mold --run \
    ninja -C ${BUILD_DIR} -v -j${THREADS} #cpp_etc dictator cpp_sim_2_1_4 phare-exe
  # (cd build && ctest -j${THREADS})
  date
) 1> >(tee $CWD/.cmake.sh.out ) 2> >(tee $CWD/.cmake.sh.err >&2 )
