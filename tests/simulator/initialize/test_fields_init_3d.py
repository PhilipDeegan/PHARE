"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
import itertools

import numpy as np
import matplotlib
from ddt import data, ddt, unpack

from pyphare.cpp import supported_particle_layouts

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 3
interp_orders = [1, 2, 3]
ppc, cells = 10, 20


def permute():
    return [
        dict(interp_order=interp_order, sim_setup_kwargs=dict(layout=layout))
        for interp_order, layout in itertools.product(
            interp_orders, supported_particle_layouts()
        )
    ]


@ddt
class Initialization3DTest(InitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, interp_order, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_B_is_as_provided_by_user(
            ndim, interp_order, nbr_part_per_cell=ppc, cells=cells, **kwargs
        )

    @data(*permute())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, interp_order, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(
            ndim, interp_order, nbr_part_per_cell=ppc, cells=cells, **kwargs
        )

    @data(*permute())
    @unpack
    def test_density_is_as_provided_by_user(self, interp_order, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_is_as_provided_by_user(
            ndim, interp_order, cells=cells, **kwargs
        )

    @data(*permute())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, interp_order, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(
            ndim, interp_order, np.asarray([20, 50, 75]), cells=10, **kwargs
        )


if __name__ == "__main__":
    unittest.main()
