#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls

set -e

FILE="tests/core/numerics/ohm/test_tileset_ohm.cpp"
FILE="tests/core/numerics/ion_updater/test_multi_updater.cpp"
CARGS=${CARGS:-""}

set -x
(
  [[ $FILE == tests/core/* ]] && (
    [ ! -f "bin/core/libphare_core.a" ] && (
        mkn clean build -Kqp test_core -d google.test,+ -a "-fPIC" ${CARGS}
        mkn clean build -Kqp core -a "-fPIC" ${CARGS} "$@"
    )
    mkn clean build -p test_core -M "${FILE}" ${CARGS} "$@" run
  ) || true

  [[ $FILE != tests/core/* ]] && (
    [ ! -f "bin/amr/libphare_amr.a" ] && (
        mkn clean build -Kqdp test_amr -a "-fPIC" ${CARGS} "$@"
    )
    mkn clean build -qp test_diagnostics -M "${FILE}" ${CARGS} "$@" run
    echo "ok"
  ) || true

) 1> >(tee $CWD/.mkn.sh.out ) 2> >(tee $CWD/.mkn.sh.err >&2 )

# stop if not CI
[ -z "$CI" ] && exit 0 # continue for more

(
  set -xe
  mkn -tp core_tests ${CARGS} clean build test run "$@"
  mkn -tp more_tests ${CARGS} clean build test run "$@"
)
