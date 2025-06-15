#!/usr/bin/env bash

set -ex
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD

CARGS=${CARGS:-""}
FILE="tests/core/numerics/ion_updater/test_updater_pp_main.cpp"
FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
TEST=" -M ${FILE} "
# CARGS="${CARGS} $(clargs)"

(
  date && rm -rf tee && mkdir tee
  for i in $(seq 44 4041); do
  (
    date && echo "SEED=$PHARE_SEED"
    PHARE_SEED=${i} mkn run -p test_core ${TEST} ${CARGS} $@ #-- --gtest_filter=IonUpdaterPPTest/2.updater
  ) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
  done
)

