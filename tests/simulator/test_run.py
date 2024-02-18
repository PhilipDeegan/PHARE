#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.run import Run

import numpy as np

import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib
from tests.simulator import SimulatorTest

cpp = cpp_lib()
startMPI()

time_step = 0.005
final_time = 1
time_step_nbr = int(final_time / time_step)
timestamps = np.arange(0, final_time+1, 0.5)
diag_dir = "phare_outputs/test_run"

def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=(100, 100),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=3,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="1",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        }
    )

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

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
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    pop = "protons"
    ph.ParticleDiagnostics(
        quantity="domain",
        write_timestamps=timestamps,
        population_name=pop,
    )
    ph.FluidDiagnostics(quantity="density", write_timestamps=timestamps, population_name=pop)

    return sim


def plot(diag_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=f"{diag_dir}/harris_divb_t{time}.png", plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )

        run.GetRanks(time).plot(
            filename=f"{diag_dir}/harris_RANKS_t{time}.png",
            plot_patches=True,
        )
        run.GetN(time, pop_name="protons").plot(
            filename=f"{diag_dir}/harris_N_t{time}.png",
            plot_patches=True,
        )
        run.GetB(time).plot(
            filename=f"{diag_dir}/harris_bx_t{time}.png",
            qty="Bx",
            plot_patches=True,
        )
        run.GetB(time).plot(
            filename=f"{diag_dir}/harris_by_t{time}.png",
            qty="By",
            plot_patches=True,
        )
        run.GetB(time).plot(
            filename=f"{diag_dir}/harris_bz_t{time}.png",
            qty="Bz",
            plot_patches=True,
        )

        run.GetJ(time).plot(
            filename=f"{diag_dir}/harris_Jz_t{time}.png",
            qty="Jz",
            plot_patches=True,
            vmin=-2,
            vmax=2,
        )



class RunTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RunTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RunTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        sim = config()
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(sim).run()
        if cpp.mpi_rank() == 0:
            plot(diag_dir)
        cpp.mpi_barrier()


if __name__ == "__main__":
    import unittest
    unittest.main()
