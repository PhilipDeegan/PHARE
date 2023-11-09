#!/usr/bin/env bash

set -ex

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD # move to project root

REPO="/mkn/r"
MPI_PATH=$(mkn -G MPI_PTH)
export LD_LIBRARY_PATH=${REPO}/google/test/master/bin/build:${MPI_PATH}:${REPO}/llnl/samrai/master/lib:${ROOT}/bin/core:${ROOT}/bin/amr:${ROOT}/bin/init
cp tests/diagnostic/job_2d.py.in bin/test_diagnostics/job_2d.py

CXXFLAGS=${CXXFLAGS:-""}
M_FILE="-M tests/diagnostic/test-diagnostics_2d.cpp"
#time (
  date
  mkn build ${M_FILE} -p test_diagnostics -a "${CXXFLAGS}" -gO 0 -x mpi $@
  (
    cd "$ROOT/bin/test_diagnostics"
    ./phare
  )
  date
#) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )
