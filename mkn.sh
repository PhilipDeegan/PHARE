#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls
set -e

export KLOG=${KLOG:="0"}
export PHARE_SCOPE_TIMING=1
PHARE_CELLS=${PHARE_CELLS:-15}
PHARE_PPC=${PHARE_PPC:-100}
TEST="-M tests/core/numerics/ion_updater/test_updater_pp_main.cpp"
ARGS="${TEST} -P mkn.base=gpu_ -x "
[ -d /opt/rocm/bin ] && ARGS+="res/mkn/hip" || ARGS+="res/mkn/clang_cuda "
set -x
( # build
  mkn clean build -p test_core ${ARGS} run $@ -- --gtest_filter=IonUpdaterPPTest/12.updater
) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
exit 0

(
  date
  rm -rf tee && mkdir tee
  export PHARE_CELLS PHARE_PPC
  for i in $(seq 1000 9999); do
    (
      PHARE_SEED=${i} mkn -M tests/core/numerics/ion_updater/test_updater_pp_main.cpp  \
          -p test_core ${ARGS} run \
          -- --gtest_filter=IonUpdaterPPTest/12.updater 2>&1 | tee "tee/test.txt" ;
      sleep 1 # no throttle

    ) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
  done
)

# which mlr &&  mlr --csv --opprint   --from   .phare_times.0.csv sort -nr dim -nr time then cut -x -f storage || echo "miller not installed"
