#!/usr/bin/env python3
"""
Tests that tile-based particle layouts produce identical initialization
and time-advance results to the reference (AoSMapped) layout.
"""

import os
import unittest
import itertools
from ddt import data, ddt, unpack

from pyphare import cpp
from pyphare.core.box import nDBox
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest


ndim_list = [1, 2, 3]
interp_orders = [1]
_ref_layout = "AoSMapped"
cells = 14
ppc_per_dim = [100, 33, 15]
# os.environ["PHARE_TILING_MIN_BEFORE_SPLIT"] = "1000"


def permute():
    cmp_layouts = [
        l for l in cpp.supported_particle_layouts() if _ref_layout not in str(l)
    ]

    return [
        dict(ndim=ndim, interp_order=interp, cmp_layout=cmp_layout)
        for ndim, interp, cmp_layout in itertools.product(
            ndim_list, interp_orders, cmp_layouts
        )
    ]


def compare_hierarchies(test, run0, run1, atol):
    """Asserts every diagnosed quantity is identical between two Runs, across
    every dumped time. `atol` is a dict with keys "b", "e", "moments", "particles".
    """
    times = run0.all_times()["B"]
    pops = run0.all_pops()

    # Priority order for debugging: particles -> moments -> B -> E, so the
    # first assertion to fail points at the earliest quantity in the causal
    # chain (particles feed moments, moments feed E, E feeds B).
    for pop in pops:
        for kind in ["levelGhost", "domain"]:
            p0 = run0.GetParticles(times, pop, type=kind)
            p1 = run1.GetParticles(times, pop, type=kind)
            test.assertTrue(hierarchy_compare(p0, p1, atol=atol["particles"]))

    Ni0 = run0.GetNi(times)
    Ni1 = run1.GetNi(times)
    test.assertTrue(hierarchy_compare(Ni0, Ni1, atol=atol["moments"]))

    Vi0 = run0.GetVi(times)
    Vi1 = run1.GetVi(times)
    test.assertTrue(hierarchy_compare(Vi0, Vi1, atol=atol["moments"]))

    for pop in pops:
        N0 = run0.GetN(times, pop_name=pop)
        N1 = run1.GetN(times, pop_name=pop)
        test.assertTrue(hierarchy_compare(N0, N1, atol=atol["moments"]))

        F0 = run0.GetFlux(times, pop_name=pop)
        F1 = run1.GetFlux(times, pop_name=pop)
        test.assertTrue(hierarchy_compare(F0, F1, atol=atol["moments"]))

    b0 = run0.GetB(times, all_primal=False)
    b1 = run1.GetB(times, all_primal=False)
    test.assertTrue(hierarchy_compare(b0, b1, atol=atol["b"]))

    e0 = run0.GetE(times, all_primal=False)
    e1 = run1.GetE(times, all_primal=False)
    test.assertTrue(hierarchy_compare(e0, e1, atol=atol["e"]))


class ALayoutInitTest(HybridInitializationTest):
    def compare_init(self, ndim, interp_order, cmp_layout, atol, **kwargs):
        common = dict(
            cells=cells,
            time_step_nbr=1,
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"init_layout_{type(self).__name__}"
        run0 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        run1 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            particle_layout=cmp_layout,
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            compare_hierarchies(self, run0, run1, atol)


@ddt
class LayoutInitL0Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, cmp_layout):
        super().compare_init(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=0, e=1e-14, moments=1e-14, particles=0),
            refinement_boxes=None,
        )

    @data(*permute())
    @unpack
    def test_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, cmp_layout)


@ddt
class LayoutInitL1Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, cmp_layout):
        super().compare_init(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=0, e=1e-14, moments=1e-14, particles=0),
            refinement_boxes={"L0": {"B0": nDBox(ndim, 1, 12)}},
        )

    @data(*permute())
    @unpack
    def test_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, cmp_layout)


class ALayoutAdvanceTest(HybridAdvanceTest):
    def compare_advance(self, ndim, interp_order, cmp_layout, atol, **kwargs):
        common = dict(
            cells=cells,
            time_step_nbr=kwargs.pop("time_step_nbr", 1),
            model_init={"seed": 1337},
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"adv_layout_{type(self).__name__}"
        run0 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        run1 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            particle_layout=cmp_layout,
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            compare_hierarchies(self, run0, run1, atol)


@ddt
class LayoutAdvanceL0Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, cmp_layout):
        super().compare_advance(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=1e-14, e=1e-12, moments=1e-14, particles=1e-14),
            refinement_boxes=None,
        )

    @data(*permute())
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, cmp_layout)


_L1_advance_atol_per_dim = {
    # deliberately at zero for now - want these to fail so the assertion
    # message's actual max diff tells us what atol they really need.
    1: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
    2: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
    3: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
}


@ddt
class LayoutAdvanceL1Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, cmp_layout):
        super().compare_advance(
            ndim,
            interp_order,
            cmp_layout,
            atol=_L1_advance_atol_per_dim[ndim],
            refinement_boxes={"L0": {"B0": nDBox(ndim, 1, 12)}},
        )

    @data(*permute())
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, cmp_layout)


# @ddt # having an L2 == regridding with refinement_boxes?
# class LayoutAdvanceL2Test(ALayoutAdvanceTest):
#     def _compare_advance(self, ndim, interp_order, cmp_layout):
#         super().compare_advance(
#             ndim,
#             interp_order,
#             cmp_layout,
#             atol=dict(b=1e-5, e=1e-5, moments=1e-5, particles=1e-5),
#             refinement_boxes={
#                 "L0": {"B0": nDBox(ndim, 1, 12)},
#                 "L1": {"B0": nDBox(ndim, 4, 21)},
#             },
#             time_step_nbr=1,
#         )

#     @data(*permute())
#     @unpack
#     def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
#         print(
#             f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
#         )
#         self._compare_advance(ndim, interp_order, cmp_layout)


if __name__ == "__main__":
    from pyphare.simulator.simulator import startMPI

    # raise "?"
    startMPI()
    unittest.main()
