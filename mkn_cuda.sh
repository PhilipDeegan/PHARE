#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls
set -ex

export KLOG=${KLOG:="0"}
PHARE_CELLS=${PHARE_CELLS:-15}
PHARE_PPC=${PHARE_PPC:-100}
ARGS=" -P mkn.base=gpu_ -x res/mkn/clang_cuda "

( # build
  mkn clean build -p test_core ${ARGS} -M tests/core/numerics/ion_updater/test_updater_pp_main.cpp run $@
) 1> >(tee $CWD/.mkn_cuda.sh.out ) 2> >(tee $CWD/.mkn_cuda.sh.err >&2 )
exit 0

(
  date
  export PHARE_SCOPE_TIMING=1
  cd "$CWD"
  rm -rf tee && mkdir tee
  export PHARE_CELLS PHARE_PPC
  for i in $(seq 1000 9999); do
    (
      PHARE_SEED=${i} mkn -M tests/core/numerics/ion_updater/test_updater_pp_main.cpp  \
          -p test_core ${ARGS} run \
          -- --gtest_filter=IonUpdaterPPTest/12.updater 2>&1 | tee "tee/test.txt" ;
      sleep 1 # no throttle

    ) 1> >(tee $CWD/.mkn_cuda.sh.out ) 2> >(tee $CWD/.mkn_cuda.sh.err >&2 )
  done
)

# which mlr &&  mlr --csv --opprint   --from   .phare_times.0.csv sort -nr dim -nr time then cut -x -f storage || echo "miller not installed"
