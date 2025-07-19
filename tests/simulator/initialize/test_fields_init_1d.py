"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
import itertools
import matplotlib
from pyphare.cpp import supported_particle_layouts
from ddt import data, ddt, unpack
from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]


def permute():
    return [
        dict(
            ndim=ndim,
            interp_order=interp_order,
            sim_setup_kwargs=dict(layout=layout),
        )
        for interp_order, layout in itertools.product(
            interp_orders, supported_particle_layouts()
        )
    ]


@ddt
class Initialization1DTest(InitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_B_is_as_provided_by_user(**kwargs)

    @data(*permute())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(**kwargs)

    @data(*permute())
    @unpack
    def test_density_is_as_provided_by_user(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_density_is_as_provided_by_user(**kwargs)

    @data(*permute())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(**kwargs)


if __name__ == "__main__":
    unittest.main()
