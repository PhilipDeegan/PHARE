#!/usr/bin/env bash
set -e
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $CWD && cd ../.. && ROOT=$PWD # move to project root
ARGS=""
EXEC="build run"
DBOPT="-gO 3"
X_FILE="" #-x clang_asan"
FILES=(
  #tests/core/data/particles/test_sorting_etc.cpp
  #tests/core/numerics/ion_updater/test_updater.cpp
  # tests/core/utilities/box/test_box.cpp
  # tests/core/utilities/sorting/test_sorting.cpp
  tests/core/numerics/pusher/test_pusher.cpp
)
time (
  date
  for FILE in ${FILES[@]}; do
    mkn $EXEC -M ${FILE} -p test_core -a "${ARGS}" $DBOPT $X_FILE $@
  done
  date
) 1> >(tee $ROOT/.mkn.sh.out ) 2> >(tee $ROOT/.mkn.sh.err >&2 )



###
#
# -r --gtest_filter=IonUpdaterTest/1.particlesUntouchedInMomentOnlyMode
###
