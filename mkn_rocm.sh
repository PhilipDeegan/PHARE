#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls
set -ex

export KLOG=${KLOG:="0"}
PHARE_CELLS=${PHARE_CELLS:-5}
PHARE_PPC=${PHARE_PPC:-888}
ARGS=" -P mkn.base=gpu_ -x res/mkn/hip_mpi -a -DPHARE_SKIP_MPI_IN_CORE"

mkn clean build -p test_core ${ARGS} -M tests/core/numerics/ion_updater/test_updater_pp_main.cpp run $@

(
  date
  # export PHARE_SCOPE_TIMING=1
  # cd "$CWD"
  # rm -rf tee && mkdir tee
  # export PHARE_CELLS PHARE_PPC
  # for i in $(seq 1000 3000); do
  #   PHARE_SEED=${i} mkn -M tests/core/numerics/ion_updater/test_updater_pp_main.cpp  \
  #       -p test_core ${ARGS} run \
  #       -- --gtest_filter=IonUpdaterPPTest/4.updater 2>&1 | tee "tee/${i}.txt" ;
  # done
)

# mlr --csv --opprint --from .phare_times.0.csv sort -nr dim -nr time then cut -x -f storage
