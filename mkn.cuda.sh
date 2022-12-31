#! /usr/bin/env bash

set -ex

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MKN_X_FILE=${MKN_X_FILE:-settings}
MKN_MPI_X_FILE=${MKN_MPI_X_FILE:-res/mkn/mpi}
MKN_GPU_X_FILE=${MKN_GPU_X_FILE:-res/mkn/clang_cuda}
export MKN_LIB_LINK_LIB=1

KLOG=2 mkn clean build -p test_core_gpu_mkn -M tests/core/data/test_vector.cpp -x ${MKN_GPU_X_FILE} -l -pthread run
