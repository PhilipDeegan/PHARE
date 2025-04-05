#! /usr/bin/env bash

set -e

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD

mkn clean build -KOp core -a -DPHARE_SKIP_MPI_IN_CORE
mkn clean build -p test_core -d google.test,+ -Og 0 -Ka -fPIC
mkn clean build -p bench -d google.benchmark,+ -Og 0 -Ka -fPIC

./tools/mkn/test.sh
./tools/mkn/bench.sh
