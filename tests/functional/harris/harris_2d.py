#!/usr/bin/env python3

from datetime import datetime

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.run import Run

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib

cpp = cpp_lib()
startMPI()

diag_outputs0 = "phare_outputs/balanceL1/harris/2d"
diag_outputs1 = "phare_outputs/balanceL1/harris/2d_rebal"

time_step = 0.001
time_step_nbr = 1
restart_time = 0
final_time = time_step * time_step_nbr

# array([10., 20., 30., 40., 50.])
timestamps = np.arange(
    restart_time, final_time + time_step, time_step * time_step_nbr / 5, dtype=float
)
timestamps = [0]


def config(diag_dig, **kwargs):
    sim = ph.Simulation(
        dl=(0.2, 0.2),
        cells=(100, 100),
        tag_buffer=1,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        # smallest_patch_size=15,
        # largest_patch_size=25,
        max_nbr_levels=2,
        refinement="tagging",
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dig, "mode": "overwrite"},
        },
        # restart_options=dict(
        #     dir=diag_dig,
        #     mode="overwrite",
        #     timestamps=timestamps,
        #     restart_time=restart_time,
        # ),
        **kwargs,
    )

    def density(x, y):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 1
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vxyz,
        "vbulky": vxyz,
        "vbulkz": vxyz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
        "nbr_part_per_cell": 100,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    pop = "protons"
    ph.ParticleDiagnostics(
        quantity="domain",
        write_timestamps=timestamps,
        population_name="protons",
    )
    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )

    ph.global_vars.sim = None
    return sim


def plot(diag_dir):
    if cpp.mpi_rank() == 0:
        run = Run(diag_dir)
        for time in timestamps:
            run.GetN(time, pop_name="protons").plot(
                filename=f"{diag_dir}/harris_N_t{time}.png",
                plot_patches=True,
            )

            for i in range(2):
                run.GetRanks(time).plot(
                    filename=f"{diag_dir}/harris_RANKS_L{i}_t{time}.png",
                    plot_patches=True,
                    levels=(i,),
                )

            run.GetRanks(time).plot(
                filename=f"{diag_dir}/harris_RANKS__t{time}.png", plot_patches=True
            )

            run.GetDivB(time).plot(
                filename=f"{diag_dir}/harris_divb_t{time}.png", plot_patches=True
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

    cpp.mpi_barrier()


def run(sim, diag_dir):
    ph.global_vars.sim = sim
    Simulator(sim).run().reset()
    plot(diag_dir)
    ph.global_vars.sim = None


def main():
    # sim0 = config(diag_outputs0)
    sim1 = config(
        diag_outputs1,
        advanced={
            "integrator/rebalance_coarsest": 1,
            "integrator/rebalance_coarsest_every": 500,
            "integrator/flexible_load_tolerance": 0.05,
        },
    )

    # run(sim0, diag_outputs0)
    run(sim1, diag_outputs1)


if __name__ == "__main__":
    main()
