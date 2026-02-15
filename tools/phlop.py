#!/usr/bin/env python3
"""
export PYTHONPATH="$HOME/git/phare/mkn:$HOME/git/phare/mkn/build:$HOME/git/phare/mkn/pyphare:$HOME/git/phlop"
./mkn.sh  -Og 0  clean build -w mkn.gpu  -a -fno-omit-frame-pointer
./tools/mkn_profile.sh
python3 tools/python3/phloping.py print_scope_timings -f .phare_times.0.txt
"""


import os
import shlex
import subprocess


def run(cmd, env=os.environ.copy()):
    subprocess.run(shlex.split(cmd), check=True, env=env)


profile = "-p test_diagnostics -x res/mkn/mpi -gO 0 -a '-O0 -fPIC' -qw mkn.gpu run  clean build"


def do(dim):
    env = os.environ.copy()
    env.update({"PHARE_LOG": "NONE", "KLOG": "0"})

    P = f"""
    mkn_r {profile} -M tests/diagnostic/test_diagnostics_{dim}d_ref.cpp
    """
    print(P.strip())
    run(P.strip(), env)

    P = f"""
    mkn_r {profile} -M tests/diagnostic/test_diagnostics_{dim}d_cmp.cpp
    """
    print(P.strip())
    run(P.strip(), env)

    P = f"""
    python3 tools/compare_diags.py phare_outputs/diags_{dim}d/ref/ phare_outputs/diags_{dim}d/cmp/
    """
    print(P.strip())
    run(P.strip(), env)

    # h5dump phare_outputs/diags_2d/cmp/ions_pop_protons_density.h5 > cmp
    # h5dump phare_outputs/diags_2d/ref/ions_pop_protons_density.h5 > ref


if __name__ == "__main__":
    do(2)
