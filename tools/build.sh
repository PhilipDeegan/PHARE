#!/usr/bin/env bash
# usage:
#  copy this file to tools/cmake.sh - and edit as you wish
#    tools/cmake.sh is ignored by git
set -ex
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root

THREADS=${THREADS:="14"}
BUILD_DIR=${BUILD_DIR:="$CWD/build"}
SAMRAI=${SAMRAI:="/mkn/r/llnl/samrai/regular"} # "" = as subproject /mkn/r/llnl/samrai/master
FFF=("${BUILD_DIR}")
CMAKE_CONFIG="-DdevMode=ON -Dasan=OFF -Dbench=OFF -DwithCaliper=OFF -DtestMPI=OFF -DtestDuringBuild=OFF -DwithKokkosTools=ON -DKokkos_ENABLE_SERIAL=ON"
CMAKE_CONFIG="${CMAKE_CONFIG} -Dphare_configurator=ON " # -DHDF5_ROOT=/usr/local/HDF_Group/HDF5/1.15.0" #  " #-DHDF5_ROOT=/opt/mpi/rocm_hdf5"
CMAKE_CONFIG="${CMAKE_CONFIG} -DwithPhlop=OFF " # -DPHARE_EXEC_LEVEL_MAX=99 "
# CMAKE_CONFIG="${CMAKE_CONFIG} -DhighResourceTests=ON " # -DPHARE_EXEC_LEVEL_MAX=99 "

# CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native"                # SUPER RELEASE
# CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer" # OPTIMZ AND DEBUG
# CMAKE_CXX_FLAGS="-g0 -O3 -march=native -mtune=native"                         # Optimz and asserts/debug
CMAKE_CXX_FLAGS="-g3 -O0 -fno-omit-frame-pointer"                             # PURE DEBUG

CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DPHARE_DIAG_DOUBLES=1 -DPHARE_LOG_LEVEL=1 " #
CMAKE_CONFIG="${CMAKE_CONFIG} -DCMAKE_BUILD_TYPE=Debug" # Or Debug RelWithDebInfo Release

exec 19>$CWD/.cmake.sh.cmd # set -x redirect
export BASH_XTRACEFD=19  # set -x redirect

CC=${CC:="gcc"}
CXX=${CXX:="g++"}
set -xe
time (
  date
  [ -n "$CLEAN" ] && (( $CLEAN == 1 )) && for f in ${FFF[@]}; do rm -rf $f; done
  [ ! -f "$CWD/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
  mkdir -p ${BUILD_DIR}
  [[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"
  (
    export CC CXX && cd ${BUILD_DIR}
    cmake $CWD ${SAMRAI} -G Ninja ${CMAKE_CONFIG} -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}"
    mold --run ninja -C ${BUILD_DIR} -v -j${THREADS} $@ # cpp cpp_etc dictator
    # ctest -j${THREADS}
    # pythom3 -m
  )
  date
) 1> >(tee $CWD/.cmake.sh.out ) 2> >(tee $CWD/.cmake.sh.err >&2 )
