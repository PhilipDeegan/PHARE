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
from pyphare.cpp import supported_particle_layouts
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

from tests.simulator.test_initialization import InitializationTest
from tests.simulator.test_advance import AdvanceTestBase


ndim_list = [1]
interp_orders = [1]
_ref_layout = "AoSMapped"
cells = 32
ppc_per_dim = [100, 66, 20]
# os.environ["PHARE_TILING_MIN_BEFORE_SPLIT"] = "1000"


def permute():
    cmp_layouts = [l for l in supported_particle_layouts() if _ref_layout not in str(l)]

    return [
        dict(ndim=ndim, interp_order=interp, cmp_layout=cmp_layout)
        for ndim, interp, cmp_layout in itertools.product(
            ndim_list, interp_orders, cmp_layouts
        )
    ]


class ALayoutInitTest(InitializationTest):
    def compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        common = dict(
            qty=qty,
            cells=cells,
            time_step_nbr=1,
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )
        ref_hier = self.getHierarchy(
            ndim,
            interp_order,
            diag_outputs=f"init_{__name__}_{_ref_layout}",
            **common,
        )
        cmp_hier = self.getHierarchy(
            ndim,
            interp_order,
            sim_setup_kwargs=dict(layout=cmp_layout),
            diag_outputs=f"init_{__name__}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            eqr = hierarchy_compare(ref_hier, cmp_hier, atol=atol)
            print(eqr)
            self.assertTrue(eqr)


@ddt
class LayoutInitL0Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0):
        super().compare_init(
            ndim, interp_order, qty, cmp_layout, atol=atol, refinement_boxes=None
        )

    @data(*permute())
    @unpack
    def test_B_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "b", cmp_layout, atol=0)

    @data(*permute())
    @unpack
    def test_E_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "e", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_moments_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "particles", cmp_layout, atol=0)


@ddt
class LayoutInitL1Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0):
        super().compare_init(
            ndim,
            interp_order,
            qty,
            cmp_layout,
            atol=atol,
            refinement_boxes={"L0": {"B0": nDBox(ndim, 5, 14)}},
        )

    @data(*permute())
    @unpack
    def test_B_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "b", cmp_layout, atol=0)

    @data(*permute())
    @unpack
    def test_E_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "e", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_moments_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "particles", cmp_layout, atol=0)


class ALayoutAdvanceTest(AdvanceTestBase):
    def compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        common = dict(
            qty=qty,
            cells=cells,
            time_step_nbr=1,
            model_init={"seed": 1337},
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"adv_layout_{type(self).__name__}"
        ref_hier = self.getHierarchy(
            ndim,
            interp_order,
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        cmp_hier = self.getHierarchy(
            ndim,
            interp_order,
            sim_setup_kwargs=dict(layout=cmp_layout),
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            eqr = hierarchy_compare(ref_hier, cmp_hier, atol=atol)
            print(eqr)
            self.assertTrue(eqr)


@ddt
class LayoutAdvanceL0Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        super().compare_advance(
            ndim, interp_order, qty, cmp_layout, atol=atol, refinement_boxes=None
        )

    @data(*permute())
    @unpack
    def test_B_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "b", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_E_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "e", cmp_layout, atol=1e-12)

    @data(*permute())
    @unpack
    def test_moments_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "particles", cmp_layout, atol=1e-14)


@ddt
class LayoutAdvanceL1Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        super().compare_advance(
            ndim,
            interp_order,
            qty,
            cmp_layout,
            atol=atol,
            refinement_boxes={"L0": {"B0": nDBox(ndim, 5, 14)}},
        )

    @data(*permute())
    @unpack
    def test_B_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "b", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_E_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "e", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_moments_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "particles", cmp_layout, atol=1e-14)


if __name__ == "__main__":
    from pyphare.simulator.simulator import startMPI

    startMPI()
    unittest.main()
