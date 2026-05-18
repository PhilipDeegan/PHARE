"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
import itertools

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box3D
from pyphare.cpp import supported_particle_layouts

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 3
interp_orders = [1, 2, 3]
ppc, cells = 10, 30


def permute(boxes={}):
    def f(interp, layout):
        dic = dict(
            ndim=ndim,
            interp_order=interp,
            sim_setup_kwargs=dict(layout=layout),
        )
        if boxes:
            return dict(refinement_boxes=boxes, **dic)
        return dic

    return [
        f(*els)
        for els in itertools.product(interp_orders, supported_particle_layouts())
    ]


@ddt
class Initialization3DTest(InitializationTest):
    @data(*permute())
    @unpack
    def test_nbr_particles_per_cell_is_as_provided(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_nbr_particles_per_cell_is_as_provided(ppc=ppc, cells=cells, **kwargs)

    @data(
        *permute({"L0": {"B0": Box3D(10, 14)}}),
        *permute({"L0": {"B0": Box3D(10, 14)}, "L1": {"B0": Box3D(22, 26)}}),
        *permute({"L0": {"B0": Box3D(2, 6), "B1": Box3D(7, 11)}}),
    )
    @unpack
    def test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, **kwargs
    ):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_levelghostparticles_have_correct_split_from_coarser_particle(
            self.getHierarchy(
                **kwargs,
                qty="particles",
                cells=cells,
                nbr_part_per_cell=ppc,
            )
        )
        print(
            f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds"
        )

    @data(
        *permute({"L0": {"B0": Box3D(10, 14)}}),
        *permute({"L0": {"B0": Box3D(5, 14)}, "L1": {"B0": Box3D(15, 19)}}),
        *permute({"L0": {"B0": Box3D(2, 12), "B1": Box3D(13, 25)}}),
    )
    @unpack
    def test_domainparticles_have_correct_split_from_coarser_particle(self, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_domainparticles_have_correct_split_from_coarser_particle(
            nbr_part_per_cell=ppc, **kwargs
        )
        print(
            f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds"
        )


if __name__ == "__main__":
    unittest.main()
