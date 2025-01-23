#!/usr/bin/env bash
# usage:
#  copy this file to tools/cmake.sh - and edit as you wish
#    tools/cmake.sh is ignored by git
set -ex
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR && cd .. && CWD=$PWD # move to project root

THREADS=${THREADS:="14"}
BUILD_DIR=${BUILD_DIR:="$CWD/build"}
SAMRAI=${SAMRAI:=""} # "" = as subproject /mkn/r/llnl/samrai/master
FFF=("${BUILD_DIR}")
# CMAKE_CONFIG="-DdevMode=ON -Dasan=OFF -Dbench=OFF -DwithCaliper=OFF -DtestMPI=OFF -DtestDuringBuild=OFF -DwithKokkosTools=ON -DKokkos_ENABLE_SERIAL=ON"
CMAKE_CONFIG="-DCMAKE_POSITION_INDEPENDENT_CODE=ON -DENABLE_OPENMP=OFF "
CMAKE_CONFIG="${CMAKE_CONFIG} -DwithKokkosTools=ON -DKokkos_ENABLE_SYCL=ON -DKokkos_ARCH_INTEL_PVC=ON -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx "
CMAKE_CONFIG="${CMAKE_CONFIG} -DMPI_GUESS_LIBRARY_NAME=OPENMPI -DMPI_EXECUTABLE_SUFFIX=.openmpi -DMPI_CXX_SKIP_MPICXX=ON -DCMAKE_CXX_COMPILER=icpx "

# CMAKE_CXX_FLAGS="-DNDEBUG -g0 -O3 -march=native -mtune=native"                # SUPER RELEASE
# CMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -fno-omit-frame-pointer" # OPTIMZ AND DEBUG
# CMAKE_CXX_FLAGS="-g0 -O3 -march=native -mtune=native"                         # Optimz and asserts/debug
CMAKE_CXX_FLAGS="-g3 -O0 -fno-omit-frame-pointer"                             # PURE DEBUG

CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DPHARE_DIAG_DOUBLES=1 -DPHARE_LOG_LEVEL=1 " #
CMAKE_CONFIG="${CMAKE_CONFIG} -DCMAKE_BUILD_TYPE=Debug" # Or Debug RelWithDebInfo Release

exec 19>$CWD/.cmake.sh.cmd # set -x redirect
export BASH_XTRACEFD=19  # set -x redirect

set -xe
time (
  date
  [ -n "$CLEAN" ] && (( $CLEAN == 1 )) && for f in ${FFF[@]}; do rm -rf $f; done
  [ ! -f "$CWD/CMakeLists.txt" ] && echo "script expected to be run from project root" && exit 1
  mkdir -p ${BUILD_DIR}
  [[ -n "${SAMRAI}" ]] && SAMRAI=-DSAMRAI_ROOT="${SAMRAI}"
  (
    cd ${BUILD_DIR}
    cmake $CWD ${SAMRAI} -G Ninja ${CMAKE_CONFIG} -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}"
    ninja -C ${BUILD_DIR} -v amper_gpu_kokkos
    ctest -R amper_gpu_kokkos -V
  )
  date
) 1> >(tee $CWD/.cmake.sh.out ) 2> >(tee $CWD/.cmake.sh.err >&2 )
