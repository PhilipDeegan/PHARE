#!/usr/bin/env python3


import pyphare.pharein as ph  # lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectronModel
import numpy as np

# configure the simulation
PPC = 75000

sim = Simulation(
    smallest_patch_size=20,
    largest_patch_size=20,
    time_step=0.001,
    time_step_nbr=1,  # number of time steps (not specified if time_step and final_time provided)
    boundary_types=["periodic"]
    * 2,  # boundary condition, string or tuple, length == len(cell) == len(dl)
    cells=[40] * 2,  # integer or tuple length == dimension
    dl=[0.3] * 2,  # mesh size of the root level, float or tuple
    # max_nbr_levels=2,  # (default=1) max nbr of levels in the AMR hierarchy
    # refinement="tagging",
    refinement_boxes={"L0": {"B0": [(15, 15), (24, 24)]}},
    diag_options={"format": "phareh5", "options": {"dir": "phare_outputs"}},
)


def density(*xyz):
    return 1.0


def bx(*xyz):
    return 1.0


def by(*xyz):
    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)


def bz(*xyz):
    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)


def vx(*xyz):
    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)


def vy(*xyz):
    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)


def vz(*xyz):
    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)


def vthx(*xyz):
    return 0.01


def vthy(*xyz):
    return 0.01


def vthz(*xyz):
    return 0.01


vvv = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthx,
    "vthy": vthy,
    "vthz": vthz,
}

MaxwellianFluidModel(
    bx=bx,
    by=by,
    bz=bz,
    protons={
        "nbr_part_per_cell":PPC,
        "charge": 1,
        "density": density,
        **vvv,
        "init": {"seed": 1337},
    },
)

ElectronModel(closure="isothermal", Te=0.12)


sim = ph.global_vars.sim

timestamps = np.arange(0, sim.final_time + sim.time_step, 100 * sim.time_step)


for quantity in ["E", "B"]:
    ph.ElectromagDiagnostics(
        quantity=quantity,
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
    )


for quantity in ["density", "bulkVelocity"]:
    ph.FluidDiagnostics(
        quantity=quantity,
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
    )

pops = [
    "protons",
]
