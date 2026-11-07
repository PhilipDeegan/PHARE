#!/usr/bin/env bash
# run by res/cmake/def.cmake at configure time (-Dphare_configurator=ON)
#  everything is generated under ${PHARE_BUILD_DIR}/phare_configurator, see config.py

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd "$CWD"
ROOT_DIR=$(cd ../.. && pwd)
[ -f "${ROOT_DIR}/res/cmake/def.cmake" ] || (echo "Corruption detected" && exit 1) # confirm correct directory

export CMAKE_COMMAND=${1}
export CMAKE_CXX_COMPILER=${2}
export PYTHON_EXECUTABLE=${3}
export PHARE_BUILD_DIR=${4}

OUT_DIR="${PHARE_BUILD_DIR}/phare_configurator"

# fresh every configure: a cached build would miss a changed compiler / MPI / HDF5
rm -rf "${OUT_DIR}" && mkdir -p "${OUT_DIR}"
"$CMAKE_COMMAND" -S "$CWD" -B "${OUT_DIR}/probe" -DCMAKE_BUILD_TYPE=Release \
                 -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}" > "${OUT_DIR}/probe.log" 2>&1 \
    || (cat "${OUT_DIR}/probe.log" && exit 1)
"$CMAKE_COMMAND" --build "${OUT_DIR}/probe" >> "${OUT_DIR}/probe.log" 2>&1 \
    || (cat "${OUT_DIR}/probe.log" && exit 1)
(cd "${OUT_DIR}" && ./probe/phare_configurator)

"$PYTHON_EXECUTABLE" config.py
