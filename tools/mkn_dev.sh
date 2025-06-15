#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls
set -e

export KLOG=${KLOG:-0}
export PHARE_SCOPE_TIMING=1

# export PHARE_PPC=${PHARE_PPC:-5}
# export PHARE_CELLS=${PHARE_CELLS:-10}
# export PHARE_PATCHES=${PHARE_PATCHES:-1}

export PHARE_ASYNC_TIMES=${PHARE_ASYNC_TIMES:-".phare/async/multi_updater"}
export PHARE_ASYNC_TIMES+="/${PHARE_PATCHES}/${PHARE_CELLS}/${PHARE_PPC}"

TEST="-M tests/core/numerics/ion_updater/test_updater_pp_main.cpp"
TEST="-M tests/core/numerics/ion_updater/test_multi_updater.cpp"
# TEST="-M tests/core/data/particles/test_particles_selecting.cpp"
ARGS="${TEST} -P mkn.base=gpu_ -x "
[ -d /opt/rocm/bin ] && ARGS+="res/mkn/hip" || ARGS+="res/mkn/clang_cuda "

set -x

# ( # build
#   [ ! -f "bin/core/libphare_core.so" ] && (
#       mkn clean build -Oqp test_core -d google.test,+ -a "-fPIC"
#       mkn clean build -tOqp core -a "-fPIC" -x res/mkn/mpi # -DPHARE_SKIP_MPI_IN_CORE
#   )
#   # mkn run -p test_core ${ARGS} $@ # -- --gtest_filter=IonUpdaterPPTest/14.updater
#   mkn run -p test_core ${ARGS} $@ # -- --gtest_filter=IonUpdaterPPTest/14.updater

#   [ -f ".phare_times.0.txt" ] && mv "${PHARE_ASYNC_TIMES}"

# )  1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
# exit 0 # comment out to do soak test

(
  date && rm -rf tee && mkdir tee
  for i in $(seq 1009 9999); do
  (
    reset
    date && echo "SEED=$PHARE_SEED"
    PHARE_SEED=${i} mkn run -p test_core ${ARGS} $@ # -- --gtest_filter=IonUpdaterPPTest/12.updater
  ) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
  done
)

# scope function times as csv
# PYTHONPATH=$PWD/subprojects/phlop python3 subprojects/phlop/tests/timing/test_scope_timer.py  test_scope_timer -f .phare_times.0.txt
# which mlr && mlr --csv --opprint --from times.csv sort -nr dim -nr time then cut -x -f storage || echo "miller not installed"

# cuda gdb example
# KLOG=3 MKN_DBG="cuda-gdb --args" PHARE_CELLS=13 PHARE_PPC=44 mkn $ARGS dbg -g -- --gtest_filter=IonUpdaterPPTest/12.updater 1> >(tee $CWD/mkn.sh.out ) 2> >(tee $CWD/mkn.sh.err >&2 )

# cuda memory checker example
# cls; MKN_DBG="compute-sanitizer --tool memcheck" PHARE_CELLS=3 PHARE_PPC=6 mkn $ARGS dbg -g -- --gtest_filter=IonUpdaterPPTest/12.updater 1> >(tee $CWD/mkn.sh.out ) 2> >(tee $CWD/mkn.sh.err >&2 )

#
# ARGS="-p test_core -M tests/core/numerics/ion_updater/test_multi_updater.cpp -P mkn.base=gpu_ -x res/mkn/clang_cuda"
# MKN_DBG="cuda-gdb" PHARE_COMPARE=0 PHARE_PPC=10 PHARE_CELLS=5 mkn $ARGS dbg
# MKN_DBG="compute-sanitizer --tool memcheck" PHARE_COMPARE=0 PHARE_PPC=10 PHARE_CELLS=5 mkn $ARGS dbg
