#!/usr/bin/env python3
import os
import sys
import numpy as np
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

np.set_printoptions(threshold=sys.maxsize)

os.environ["PHARE_SCOPE_TIMING"] = "0"  # turn on scope timing


ph.NO_GUI()


cells = (101, 101, 101)
# cells = (41, 41, 41)
dl = (0.1, 0.1, 0.1)
dx, dy, dz = dl

name = "bowler"
diag_outputs = f"phare_outputs/test/{name}"

final_time = 35
time_step = 0.001
timestamps = [final_time]


def b3(sim, x, y, z):
    L = sim.simulation_domain()[0]
    mid = L / 2

    X = x - mid
    Y = y - mid
    Z = z - mid

    # U, V, W = -X, -Y, -Z

    # # Normalize vectors to unit length
    # magnitude = np.sqrt(U**2 + V**2 + W**2) + 1e-5  # Avoid division by zero

    # # Normalize vectors (unit length)
    # U /= magnitude
    # V /= magnitude
    # W /= magnitude

    # # Define circular mask (radius = 25)
    # radius = L * 0.4
    # diff = 0.2

    # # # outer mask
    # # mask = X**2 + Y**2 + Z**2 <= (radius + diff) ** 2
    # # U[~mask] = 0
    # # V[~mask] = 0
    # # W[~mask] = 0

    # # # inner mask
    # # mask = X**2 + Y**2 + Z**2 >= (radius - diff) ** 2
    # # U[~mask] = 0
    # # V[~mask] = 0
    # # W[~mask] = 0

    # # anything outside the horizontal middle plane = 0
    # # mask = np.abs(Y) < 5  # & (X**2 + Z**2 >= (radius - diff) ** 2)
    # # mask = np.abs(Y) < 5 & (X**2 + Z**2 <= radius**2)

    # mask = (np.abs(Y) < dy / 2) & (
    #     (np.abs(np.abs(X) - mid) < dx) | (np.abs(np.abs(Z) - mid) < dz)
    # )

    # U[~mask] = 0
    # V[~mask] = 0
    # W[~mask] = 0

    wavelength1 = 2.0

    k1 = 2 * np.pi / wavelength1

    frequency1 = 1.0

    omega1 = 2 * np.pi * frequency1

    # Radial distance
    R = np.sqrt(X**2 + Y**2 + Z**2)

    # Radial inward direction (unit vectors)
    eps = 1e-8
    Ux_dir = -X / (R + eps)
    Uy_dir = -Y / (R + eps)
    Uz_dir = -Z / (R + eps)

    # Two inward traveling radial waves
    t = 0
    phi1 = k1 * R + omega1 * t

    A = np.cos(phi1)

    # Final vector field
    U = A * Ux_dir
    V = A * Uy_dir
    W = A * Uy_dir

    U *= 0.001
    V *= 0.001
    W *= 0.001

    return U, V, W


_globals = dict(ts=0)


def update(postOp):
    from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim

    live = postOp.live
    sim = live.simulation

    _globals["ts"] += 1
    ts = _globals["ts"]

    # print("ts", ts)
    if ts % 100 != 0:
        return
    # print("ts++", ts)

    hier = None
    for i, c in enumerate(["x", "y", "z"]):
        hier = hierarchy_from_sim(live, qty=f"EM_B_{c}", hier=hier)
    # for lvl_nbr, level in hier.levels(hier.times()[0]).items():
    L0 = hier.level(0, hier.times()[0])
    for ip, patch in enumerate(L0.patches):
        pdata_names = list(patch.patch_datas.keys())
        for i, name in enumerate(pdata_names):
            pd = patch.patch_datas[name]
            nbrGhosts = pd.ghosts_nbr
            select = tuple([slice(nbrGhost, -(nbrGhost)) for nbrGhost in nbrGhosts])
            pd[pd.box] += b3(sim, *pd.meshgrid())[i][select]


def config():
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        dl=dl,
        cells=cells,
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "pharevtkhdf",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        restart_options={
            "dir": "checkpoints",
            "mode": "overwrite",
            # "elapsed_timestamps": [0],
            "timestamps": [15, 20, 25, 30, final_time],
            "restart_time": "auto",
        },
        strict="very",
    )

    def density(x, y, z):
        return 0.5

    def b(x, y, z):
        return b3(sim, x, y, z)

    def bx(x, y, z):
        return b(x, y, z)[0]

    def by(x, y, z):
        return b(x, y, z)[1]

    def bz(x, y, z):
        return b(x, y, z)[2]

    def vxyz(x, y, z):
        return 0.0

    def vthxyz(x, y, z):
        return 0.00001

    C = "xyz"
    vvv = {
        **{f"vbulk{c}": vxyz for c in C},
        **{f"vth{c}": vthxyz for c in C},
        "nbr_part_per_cell": 100,
    }
    protons = {
        "charge": 1,
        "mass": 1,
        "density": density,
        **vvv,
        "init": {"seed": 12334},
    }
    ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["mass_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    # ph.InfoDiagnostics(quantity="particle_count")  # defaults all coarse time steps

    return sim


if __name__ == "__main__":
    Simulator(config()).run()
    main()
