#!/usr/bin/env bash
CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && cd $CWD/.. && CWD=$PWD
set -o pipefail
shopt -s expand_aliases
alias cls="clear; printf '\033[3J'"
cls
set -e

export PHARE_GPU_BYTES=5000000000 # 5gb
export PHARE_SCOPE_TIMING=1
export PHARE_ASYNC_THREADS=3

TEST="-M tests/core/numerics/ion_updater/test_multi_updater.cpp"
ARGS="${TEST} -P mkn.base=gpu_ -x "
[ -d /opt/rocm/bin ] && ARGS+="res/mkn/hip" || ARGS+="res/mkn/clang_cuda "

(
    python3 -O tools/mkn_profile.py $ARGS $@
)
