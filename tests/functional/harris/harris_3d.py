#!/usr/bin/env python3

import numpy as np

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib

from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI


ph.NO_GUI()
cpp = cpp_lib()

start_time = 0.0

cells = (100, 50, 25)
dl = (0.4, 0.4, 0.4)
diag_outputs = "phare_outputs/harris_3d"
time_step = 0.002
final_time = 0.002  # 10
timestamps = [0, final_time]

hs = hour_seconds = 3600.0
elapsed_restart_timestamps = [hs * 3, hs * 6, hs * 9, hs * 12, hs * 15]
ppc = 10


def diag_timestamps():
    dt = 1000 * time_step
    nt = (final_time - start_time) / dt
    return start_time + dt * np.arange(nt)


timestamps = diag_timestamps()


def config():
    L = 0.5

    sim = ph.Simulation(
        # dry_run=1,
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_nbr_levels=3,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="10",
        tagging_threshold=0.4,
        hyper_resistivity=0.008,
        hyper_mode="constant",
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": ".", "mode": "overwrite"},
        },
        # restart_options={
        #     "dir": "checkpoints",
        #     "mode": "overwrite",
        #     "elapsed_timestamps": elapsed_restart_timestamps,
        #     # "restart_time":start_time
        # },
    )

    def density(x, y, z):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y, z):
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

    def bx(x, y, z):
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

    def bz(x, y, z):
        return 0.0

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    def T(x, y, z):
        K = 0.7
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y, z):
        return 0.0

    def vy(x, y, z):
        return 0.0

    def vz(x, y, z):
        return 0.0

    def vthx(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthy(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthz(x, y, z):
        return np.sqrt(T(x, y, z))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": ppc,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "init": {"seed": cpp.mpi_rank() + 12},
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.LoadBalancer(active=True, mode="nppc", tol=0.05, every=1000)

    pop = "protons"
    ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=timestamps)
    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name=pop
    )

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    # timestamps = 23.0 + np.arange(100) * sim.time_step
    # for name in ["domain", "levelGhost"]:
    #     ph.ParticleDiagnostics(
    #         quantity=name, write_timestamps=timestamps, population_name="protons"
    #     )

    ph.InfoDiagnostics(quantity="particle_count", write_timestamps=timestamps)

    return sim


if __name__ == "__main__":
    startMPI()
    Simulator(config(), print_one_line=False).run()
