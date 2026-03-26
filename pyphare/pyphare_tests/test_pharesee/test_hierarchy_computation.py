import unittest
import numpy as np
from ddt import ddt

from pyphare.core.box import Box2D
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.core.operators import dot, cross, sqrt, modulus, grad

diag_outputs = "phare_outputs/test_hierarchy_computation"
timestamps = [0]


@ddt
class PatchHierarchyTest(unittest.TestCase):
    def setUp(self):
        import pyphare.pharein as ph

        def config():
            sim = ph.Simulation(
                final_time=0.001,
                time_step=0.001,
                cells=(50, 50),
                dl=(0.1, 0.1),
                refinement_boxes={"L0": [Box2D(10, 39)], "L1": [Box2D(30, 49)]},
                hyper_resistivity=0.005,
                resistivity=0.001,
                diag_options={
                    "format": "phareh5",
                    "options": {"dir": diag_outputs, "mode": "overwrite"},
                },
            )

            def density(x, y):
                L = sim.simulation_domain()[1]
                return 0.1

            def bx(x, y):
                return 0.1

            def by(x, y):
                return 0

            def bz(x, y):
                return 0.0

            def vx(x, y):
                return 0.0

            def vy(x, y):
                return 0.0

            def vz(x, y):
                return 0.0

            def vthx(x, y):
                return 1e-10

            def vthy(x, y):
                return 1e-10

            def vthz(x, y):
                return 1e-10

            vvv = {
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbr_part_per_cell": 100,
            }

            ph.MaxwellianFluidModel(
                bx=bx,
                by=by,
                bz=bz,
                protons={
                    "charge": 1,
                    "density": density,
                    **vvv,
                    "init": {"seed": 12334},
                },
            )

            ph.ElectronModel(closure="isothermal", Te=0.1)

            for quantity in ["E", "B"]:
                ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

            for quantity in ["charge_density", "mass_density"]:
                ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

            return sim

        Simulator(config()).initialize().reset()

    def _test_patch_hierarchy_selection_interpolation(self):
        run = Run(diag_outputs)
        Ni = run.GetNi(time=0.0)
        domain = Ni[:]
        print(domain)

    def test_all(self):  # DO NOT RUN MULTIPLE SIMULATIONS!
        checks = 0
        for test in [method for method in dir(self) if method.startswith("_test_")]:
            getattr(self, test)()
            checks += 1
        self.assertEqual(checks, 1)  # update if you add new tests


if __name__ == "__main__":
    unittest.main()
