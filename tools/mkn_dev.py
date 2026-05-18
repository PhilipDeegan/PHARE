"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import os
import shlex
import subprocess
import unittest
import itertools
import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box1D
from pyphare.cpp import supported_particle_layouts

matplotlib.use("Agg")  # for systems without GUI

ppc = 5
cells = 14
test = "-M tests/core/numerics/ion_updater/test_multi_updater.cpp"


def run(cmd, env=os.environ.copy()):
    subprocess.run(shlex.split(cmd), check=True, env=env)


def do(seed):
    env = os.environ.copy()
    env.update(
        {
            "PHARE_LOG": "NONE",
            # "PHARE_CELLS": f"{cells}",
            # "PHARE_PPC": f"{ppc}",
            "KLOG": "0",
            "PHARE_SEED": f"{seed}",
        }
    )

    P = f"mkn {test} run -p test_core"
    run(P.strip(), env)


def permute():
    return [dict(seed=i) for i in range(9999)]


@ddt
class SoakTest(unittest.TestCase):
    @data(*permute())
    @unpack
    def test_many_seeds(self, seed):
        do(seed)


if __name__ == "__main__":
    unittest.main()
