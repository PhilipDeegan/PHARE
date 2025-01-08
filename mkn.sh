#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

# TEST="-M tests/core/numerics/ion_updater/test_multi_updater.cpp"
# TEST="-M tests/core/data/tiles/test_tile.cpp"
TEST="-M tests/core/data/particles/test_particles_selecting.cpp"

XFILE="${XFILE:-res/mkn/clang_cuda}"
[ -d /opt/rocm/bin ] && XFILE="res/mkn/hip"

ARGS="${TEST} -P mkn.base=gpu_"
[ -n "XFILE" ] && ARGS+=" -x ${XFILE}"

set -x

(
  [ ! -f "bin/core/libphare_core.so" ] && (
      mkn clean build -Oqgp test_core -d google.test,+ -a "-fPIC"
      mkn clean build -tOqgp core -a "-fPIC" -x res/mkn/mpi
  )

  mkn -p test_core ${ARGS} $@

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )
