#
# parsing PHARE scope funtion timers
#

import os
import shlex
import argparse
import numpy as np
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import hierarchy_from


def subprocess_run(cmd, on_error="error message"):
    try:
        proc = subprocess.run(shlex.split(cmd), check=True, capture_output=True)
        return (proc.stdout + proc.stderr).decode().strip()
    except Exception:
        return on_error


def mkn_rocm_available():
    V = 0
    # [ -d /opt/rocm/bin ] && V=$((V + 1))
    # (( V == 0 )) && which hipcc 2>&1 > /dev/null && V=1
    # echo $V


def mkn_cuda_available():
    V = 0
    # which clang 2>&1 > /dev/null && clang -v 2>&1 | grep -q "Found CUDA" && V=1
    # echo $V


def print_env_settings():
    print("lol")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("usage: $function_name -h")
        print(
            "available functions:\n\t"
            + "\n\t".join([k for k, v in globals().items() if k.startswith("print_")]),
        )
    elif len(sys.argv) > 1:
        fn = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        if fn not in globals():
            raise ValueError("requested function does not exist")
        globals()[fn]()
