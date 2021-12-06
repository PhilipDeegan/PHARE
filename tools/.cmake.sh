#!/usr/bin/env bash

# Copy this file to $PWD/cmake.sh - this is ignored via git so you can edit it however you like
#

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && ROOT=$PWD # move to project root
exec 19>$ROOT/.cmake.sh.cmd # set -x redirect
export BASH_XTRACEFD=19  # set -x redirect
THREADS=${THREADS:="$(nproc --all)"}
BUILD_DIR=${BUILD_DIR:="$ROOT/build"}
SAMRAI=${SAMRAI:=""} # "" = check system and if not found as subproject
FFF=("${BUILD_DIR}")
NINJA=${NINJA:=""} # additional ninja flags

CMAKE_OPTS="-G Ninja -DdevMode=ON -Dasan=OFF -Dbench=OFF -DPHARE_EXEC_LEVEL_MAX=8"
CMAKE_CXX_FLAGS="-g3 -O3 -DPHARE_DIAG_DOUBLES=1 -Wfatal-errors" # -march=native -mtune=native
[[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"
echo cmake $ROOT ${SAMRAI} ${CMAKE_OPTS} -DCMAKE_CXX_FLAGS=\"${CMAKE_CXX_FLAGS}\" > .cmake.cmd

set -xe
time (
  date
  # export CC=clang-11 CXX=clang++-11
  [ -n "$CLEAN" ] && (( $CLEAN == 1 )) && for f in ${FFF[@]}; do rm -rf $f; done
  [ ! -f "$ROOT/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
  [ ! -d "$ROOT/subprojects/cppdict/include" ] && git submodule update --init
  mkdir -p ${BUILD_DIR}
  (
    cd ${BUILD_DIR}
    [[ -f "$BUILD_DIR/CMakeCache.txt" ]] && echo "cmake already configured" || \
        cmake $ROOT ${SAMRAI} ${CMAKE_OPTS} -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}"
  )
  ninja -C ${BUILD_DIR} -k 1 -v -j${THREADS} ${NINJA} || exit 1 #cpp_etc dictator cpp_sim_2_1_4
  # (cd ${BUILD_DIR} && ctest -j11)
  date
) 1> >(tee $ROOT/.cmake.sh.out ) 2> >(tee $ROOT/.cmake.sh.err >&2 )
