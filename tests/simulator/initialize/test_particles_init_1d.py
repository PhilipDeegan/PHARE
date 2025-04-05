"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box1D
from pyphare.cpp import cpp_lib
from pyphare.cpp import supported_particle_layouts

from tests.simulator.test_initialization import InitializationTest

cpp = cpp_lib()

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]
layouts = [1, 3]
ppc = 25


def permute(boxes={}):
    import itertools

    def f(interp, layout):
        dic = dict(ndim=ndim, interp_order=interp, sim_setup_kwargs=dict(layout=layout))
        if boxes:
            return dict(refinement_boxes=boxes, **dic)
        return dic

    return [
        f(*els)
        for els in itertools.product(interp_orders, supported_particle_layouts())
    ]


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class Initialization1DTest(InitializationTest):
    @data(*permute())
    @unpack
    def test_nbr_particles_per_cell_is_as_provided(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_nbr_particles_per_cell_is_as_provided(**kwargs)

    @data(
        *permute(({"L0": {"B0": Box1D(10, 14)}})),
        *permute(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 35)}})),
        *permute(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
    )
    @unpack
    def test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, **kwargs
    ):
        print(f"{self._testMethodName}_{ndim}d")

        self._test_levelghostparticles_have_correct_split_from_coarser_particle(
            self.getHierarchy(**kwargs, qty="particles", cells=30)
        )

    @data(
        *permute(({"L0": {"B0": Box1D(10, 14)}})),
        *permute(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 35)}})),
        *permute(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
    )
    @unpack
    def test_domainparticles_have_correct_split_from_coarser_particle(self, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")

        self._test_domainparticles_have_correct_split_from_coarser_particle(**kwargs)

    @data({"cells": 40, "smallest_patch_size": 20, "largest_patch_size": 20})
    def test_no_patch_ghost_on_refined_level_case(self, simInput):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_patch_ghost_on_refined_level_case(ndim, False, **simInput)

    @data({"cells": 40, "interp_order": 1})
    def test_has_patch_ghost_on_refined_level_case(self, simInput):
        print(f"{self._testMethodName}_{ndim}d")
        from pyphare.pharein.simulation import check_patch_size

        _, smallest_patch_size = check_patch_size(ndim, **simInput)
        simInput["smallest_patch_size"] = smallest_patch_size
        simInput["largest_patch_size"] = smallest_patch_size
        self._test_patch_ghost_on_refined_level_case(ndim, True, **simInput)

    @data("berger", "tile")
    def test_amr_clustering(self, clustering):
        dim = 1
        interp_order = 1
        self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": [(10,), (20,)]}},
            "particles",
            clustering=clustering,
        )


if __name__ == "__main__":
    unittest.main()
