import copy

import time
import datetime
import unittest
import numpy as np
from pathlib import Path
from datetime import timedelta

from ddt import ddt, data, unpack

from pyphare.cpp import cpp_lib

cpp = cpp_lib()

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags
from pyphare.pharesee.hierarchy.patchdata import ParticleData
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5

from pybindlibs import dictator


def permute(dic, expected_num_levels):
    # from pyphare.pharein.simulation import supported_dimensions # eventually
    dims = [1]  # supported_dimensions()
    return [
        [dim, interp, dic, expected_num_levels] for dim in dims for interp in [1, 2, 3]
    ]


def setup_model(ppc=100):
    def density(x):
        return 1.0

    def bx(x):
        return 0.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(x):
        L = ph.global_vars.sim.simulation_domain()[0]
        v1, v2 = -1, 1.0
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    def bz(x):
        return 0.5

    def b2(x):
        return bx(x) ** 2 + by(x) ** 2 + bz(x) ** 2

    def T(x):
        K = 1
        return 1 / density(x) * (K - b2(x) * 0.5)

    def vx(x):
        return 2.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vxalpha(x):
        return 3.0

    def vthxyz(x):
        return T(x)

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
    }
    vvvalpha = {
        "vbulkx": vxalpha,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
    }
    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "mass": 1,
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
        alpha={
            "mass": 4.0,
            "charge": 1,
            "density": density,
            **vvvalpha,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model


timestep = 0.001
out = "phare_outputs/restarts"
simArgs = dict(
    time_step_nbr=1,
    time_step=timestep,
    boundary_types="periodic",
    cells=200,
    dl=0.3,
    diag_options=dict(format="phareh5", options=dict(dir=out, mode="overwrite")),
    restart_options=dict(dir=out, mode="overwrite"),
)


def dup(dic={}):
    dic.update(copy.deepcopy(simArgs))
    return dic


@ddt
class RestartsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RestartsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RestartsTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def check_diags(self, diag_dir0, diag_dir1, pops, timestamps, expected_num_levels):
        if cpp.mpi_rank() > 0:
            return

        def count_levels_and_patches(qty):
            n_levels = len(qty.levels())
            n_patches = 0
            for ilvl in qty.levels().keys():
                n_patches += len(qty.level(ilvl).patches)
            return n_levels, n_patches

        self.assertGreater(len(timestamps), 0)
        for time in timestamps:
            checks = 0

            run0 = Run(diag_dir0)

            datahier0 = get_all_available_quantities_from_h5(diag_dir0, time)
            datahier1 = get_all_available_quantities_from_h5(diag_dir1, time)

            self.assertEqual(
                set(datahier0.quantities()),
                set(datahier1.quantities()),
            )

            self.assertEqual(len(datahier0.levels()), len(datahier1.levels()))
            for ilvl in range(len(datahier0.levels())):
                self.assertEqual(
                    len(datahier0.level(ilvl).patches),
                    len(datahier1.level(ilvl).patches),
                )
                for patch0, patch1 in zip(
                    datahier0.level(ilvl).patches, datahier1.level(ilvl).patches
                ):
                    self.assertEqual(patch0.box, patch1.box)

            self.assertGreater(len(datahier0.levels()), 0)

            for ilvl, lvl0 in datahier0.levels().items():
                patch_level1 = datahier1.levels()[ilvl]
                for p_idx, patch0 in enumerate(lvl0):
                    patch1 = patch_level1.patches[p_idx]
                    for pd_key, pd0 in patch0.patch_datas.items():
                        pd1 = patch1.patch_datas[pd_key]
                        self.assertNotEqual(id(pd0), id(pd1))
                        self.assertEqual(type(pd0), type(pd1))
                        if isinstance(pd1, ParticleData):
                            try:
                                self.assertEqual(pd0.dataset, pd1.dataset)
                            except AssertionError:
                                print(
                                    f"FAILED domain particles at time {time} {ilvl} {patch1.box} {patch0.box}"
                                )
                        else:
                            np.testing.assert_equal(pd0.dataset[:], pd1.dataset[:])
                        checks += 1

            n_levels, n_patches = count_levels_and_patches(
                run0.GetB(time, all_primal=False)
            )
            self.assertEqual(n_levels, expected_num_levels)
            self.assertGreaterEqual(n_patches, n_levels)  # at least one patch per level

    @data(
        *permute(
            dup(
                dict(
                    max_nbr_levels=3,
                    refinement="tagging",
                )
            ),
            expected_num_levels=3,
        ),
        *permute(dup(dict()), expected_num_levels=2),  # refinement boxes set later
    )
    @unpack
    def test_restarts(self, ndim, interp, simInput, expected_num_levels):
        print(f"test_restarts dim/interp:{ndim}/{interp}")

        simput = copy.deepcopy(simInput)

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * ndim

        if "refinement" not in simput:
            # three levels has issues with refinementboxes and possibly regridding
            b0 = [[10] * ndim, [19] * ndim]
            simput["refinement_boxes"] = {"L0": {"B0": b0}}

        # if restart time exists it "loads" from restart file
        #  otherwise just saves restart files based on timestamps
        assert "restart_time" not in simput["restart_options"]

        simput["interp_order"] = interp
        time_step = simput["time_step"]
        time_step_nbr = simput["time_step_nbr"]

        restart_idx = 0
        restart_time = time_step * restart_idx
        timestamps = [restart_time, time_step * time_step_nbr]

        # first simulation
        local_out = self.unique_diag_dir_for_test_case(f"{out}/test", ndim, interp)
        simput["restart_options"]["dir"] = local_out
        simput["restart_options"]["timestamps"] = [restart_time]
        simput["diag_options"]["options"]["dir"] = local_out
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" not in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        # self.register_diag_dir_for_cleanup(local_out)
        diag_dir0 = local_out

        # second restarted simulation
        local_out = f"{local_out}_n2"
        simput["diag_options"]["options"]["dir"] = local_out
        simput["restart_options"]["restart_time"] = restart_time
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        sim = Simulator(ph.global_vars.sim).setup()
        dictator.dump_dict(f"{local_out}/dict.data")
        sim.reset()
        # self.register_diag_dir_for_cleanup(local_out)
        # diag_dir1 = local_out

        # self.check_diags(
        #     diag_dir0, diag_dir1, model.populations, timestamps, expected_num_levels
        # )


if __name__ == "__main__":
    unittest.main()
