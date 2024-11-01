#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

TEST="-M tests/core/numerics/ion_updater/test_multi_updater.cpp"
ARGS="${TEST} -P mkn.base=gpu_ -x "
[ -d /opt/rocm/bin ] && ARGS+="res/mkn/hip" || ARGS+="res/mkn/clang_cuda "

set -x

(
  [ ! -f "bin/core/libphare_core.so" ] && (
      mkn clean build -Oqp test_core -d google.test,+ -a "-fPIC"
      mkn clean build -tOqp core -a "-fPIC" -x res/mkn/mpi
  )

  mkn clean build -p test_core ${ARGS} $@

) #1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

exit 0 # comment out to do soak test

