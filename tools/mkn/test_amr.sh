#!/usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD # move to project root

CXXFLAGS=${CXXFLAGS:-""}
M_FILE="-M tests/amr/data/particles/copy_overlap/test_particledata_copy_periodicNd.cpp"
#time (
  date
  mkn build dbg ${M_FILE} -p test_amr -a "${CXXFLAGS}" -gO 0 -x mpi $@
  date
#) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )
