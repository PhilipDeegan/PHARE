"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
import itertools
import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box1D
from pyphare.cpp import supported_particle_layouts

from tests.simulator.test_advance import AdvanceTestBase

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]


def permute(boxes={}):
    return [
        dict(
            interp_order=interp_order,
            refinement_boxes=boxes,
            sim_setup_kwargs=dict(layout=layout),
        )
        for interp_order, layout in itertools.product(
            interp_orders, supported_particle_layouts()
        )
    ]


@ddt
class Advance1DTest(AdvanceTestBase):
    @data(
        *permute({}),
        *permute({"L0": [Box1D(10, 19)]}),
        *permute({"L0": [Box1D(8, 20)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal(self, interp_order, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr = 3
        time_step = 0.001

        datahier = self.getHierarchy(
            ndim,
            interp_order,
            qty="eb",
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            **kwargs,
        )
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)

    @data(
        *permute({}),
        *permute({"L0": [Box1D(10, 19)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(
        self, interp_order, **kwargs
    ):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr = 3
        time_step = 0.001
        from pyphare.pharein.simulation import check_patch_size

        largest_patch_size, smallest_patch_size = check_patch_size(
            ndim, interp_order=interp_order, cells=[60] * ndim
        )
        datahier = self.getHierarchy(
            ndim,
            interp_order,
            qty="eb",
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=smallest_patch_size,
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            **kwargs,
        )
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)

    @data(
        *permute(({"L0": {"B0": Box1D(10, 19)}})),
        *permute(({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}})),
        *permute(({"L0": {"B0": Box1D(6, 23)}})),
        *permute(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
        *permute(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 19)}})),
        *permute(
            (
                {
                    "L0": {"B0": Box1D(5, 20)},
                    "L1": {"B0": Box1D(12, 38)},
                    "L2": {"B0": Box1D(30, 52)},
                }
            )
        ),
    )
    @unpack
    def test_field_coarsening_via_subcycles(self, interp_order, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_field_coarsening_via_subcycles(ndim, interp_order, **kwargs)

    @unittest.skip("should change to work on moments")
    @data(  # only supports a hierarchy with 2 levels
        *permute(({"L0": [Box1D(5, 9)]})),
        *permute(({"L0": [Box1D(5, 24)]})),
        *permute(({"L0": [Box1D(5, 9), Box1D(20, 24)]})),
    )
    @unpack
    def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, interp_order, **kwargs
    ):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
            ndim, interp_order, **kwargs
        )


if __name__ == "__main__":
    unittest.main()
