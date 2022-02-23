



from pyphare.cpp import cpp_lib
cpp = cpp_lib()

from tests.diagnostic import dump_all_diags

from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy import hierarchy_from, h5_filename_from, h5_time_grp_key
import pyphare.pharein as ph
import unittest
import os
import h5py
import numpy as np
from ddt import ddt, data
from pyphare.pharesee.run import Run


from tests.simulator import SimulatorTest


def setup_model(ppc=100):

    def density(x):
        return 1.

    def S(x,x0,l):
        return 0.5*(1+np.tanh((x-x0)/l))

    def bx(x):
        return 0.

    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[0]
        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))

    def bz(x):
        return 0.5

    def b2(x):
        return bx(x)**2 + by(x)**2 + bz(x)**2

    def T(x):
        K = 1
        return 1/density(x)*(K - b2(x)*0.5)

    def vx(x):
        return 2.

    def vy(x):
        return 0.

    def vz(x):
        return 0.

    def vthx(x):
        return T(x)

    def vthy(x):
        return T(x)

    def vthz(x):
        return T(x)

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    model = ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"mass":1, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 1337}},
        alpha={"mass":4, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 2334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/restarts_test/"
simArgs = dict(
  max_nbr_levels=2,
  refinement="tagging",
  time_step_nbr = 6,
  time_step = 0.001,
  boundary_types = "periodic",
  cells = 200,
  dl = 0.3,
  diag_options = {"format": "phareh5", "options": {"dir": out}},
  restart_options = {"format": "phareh5", "options": {"dir": out}},
)

def dup(dic):
    dic.update(simArgs.copy())
    return dic


@ddt
class RestartsTest(SimulatorTest):

    _test_cases = (
      dup(dict(
        smallest_patch_size=10,
        largest_patch_size =20,)),
      dup(dict(
        smallest_patch_size=20,
        largest_patch_size =20,)),
      dup(dict(
        smallest_patch_size=20,
        largest_patch_size =40,)),
    )

    def __init__(self, *args, **kwargs):
        super(RestartsTest, self).__init__(*args, **kwargs)
        self.simulator = None





    def tearDown(self):
        super(RestartsTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None


    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]


    @data(*_test_cases)
    def test_dump_diags(self, simInput):
        for ndim in [1]: #supported_dimensions(): # change when 2d tagging is in
            self._test_dump_diags(ndim, **simInput)

    def _test_dump_diags(self, dim, **simInput):
        test_id = self.ddt_test_id()

        time_step = simInput["time_step"]
        time_step_nbr = simInput["time_step_nbr"]

        # configure simulation dim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(dim)]

        if "refinement" not in simInput:
            b0 = [[10 for i in range(dim)], [19 for i in range(dim)]]
            simInput["refinement_boxes"] = {"L0": {"B0": b0}}


        py_attrs = [f"{dep}_version" for dep in ["samrai", "highfive", "pybind"] ]
        py_attrs += ["git_hash"]

        for interp in [1]: #,2,3]:
            print("test_dump_diags dim/interp:{}/{}".format(dim, interp))

            local_out = f"{out}_dim{dim}_interp{interp}_mpi_n_{cpp.mpi_size()}_id{test_id}"

            diag_dir0 = local_out

            simInput["restart_options"]["options"]["dir"] = local_out
            simInput["diag_options"]["options"]["dir"] = local_out
            self.register_diag_dir_for_cleanup(local_out)

            simulation = ph.Simulation(**simInput)
            self.assertTrue(len(simulation.cells) == dim)
            model = setup_model()

            restart_idx = 4
            restart_time=time_step * restart_idx
            ph.Restarts(write_timestamps=np.array([restart_time]))

            timestamps = [time_step * restart_idx, time_step * time_step_nbr]
            dump_all_diags(model.populations, timestamps=np.array(timestamps))

            Simulator(simulation).run().reset()
            ph.global_vars.sim = None


            local_out = f"{local_out}_n2"
            diag_dir1 = local_out
            simInput["diag_options"]["options"]["dir"] = local_out
            self.register_diag_dir_for_cleanup(local_out)

            simulation = ph.Simulation(**simInput)
            model = setup_model()
            dump_all_diags(model.populations, timestamps=np.array(timestamps))

            simulation.restart_options["options"]["restart_time"] = restart_time
            simulation.restart_options["options"]["restart_idx"] = restart_idx
            Simulator(simulation).run(restart_time=restart_time).reset()


            def check(qty0, qty1, checker):
                checks = 0
                for ilvl, lvl0 in qty0.patch_levels.items():
                    patch_level1 = qty1.patch_levels[ilvl]
                    for p_idx, patch0 in enumerate(lvl0):
                        patch1 = patch_level1.patches[p_idx]
                        for pd_key, pd0 in patch0.patch_datas.items():
                            pd1 = patch1.patch_datas[pd_key]
                            self.assertTrue(id(pd0) != id(pd1))
                            checker(pd0, pd1)
                            checks += 1
                return checks

            def check_particles(qty0, qty1):
                return check(qty0, qty1, lambda pd0, pd1: self.assertEqual(pd0.dataset, pd1.dataset))

            def check_field(qty0, qty1):
                return  check(qty0, qty1, lambda pd0, pd1: np.testing.assert_equal(pd0.dataset[:], pd1.dataset[:]))


            pops = [*model.populations]
            for time in timestamps:
                checks = 0

                run0 = Run(diag_dir0)
                run1 = Run(diag_dir1)
                checks += check_particles(run0.GetParticles(time, pops), run1.GetParticles(time, pops))
                checks += check_field(run0.GetB(time), run1.GetB(time))
                checks += check_field(run0.GetE(time), run1.GetE(time))
                checks += check_field(run0.GetNi(time), run1.GetNi(time))
                checks += check_field(run0.GetVi(time), run1.GetVi(time))

                for pop in pops:
                    checks += check_field(run0.GetFlux(time, pop), run1.GetFlux(time, pop))
                    checks += check_field(run0.GetN(time, pop), run1.GetN(time, pop))

                self.assertTrue(checks >= 14)

            ph.global_vars.sim = None



if __name__ == "__main__":
    unittest.main()


