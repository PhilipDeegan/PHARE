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

ppc = 100
cells = (101, 101)
# cells = (41, 41, 41)
dl = (0.1, 0.1)
dx, dy = dl

name = "bowler"
diag_outputs = f"phare_outputs/test/{name}"

final_time = 122
time_step = 0.001
timestamps = [final_time]


def b2(sim, x, y):
    L = sim.simulation_domain()[0]
    mid = L / 2

    X = x - mid
    Y = y - mid

    wavelength1 = 3

    k1 = 2 * np.pi / wavelength1

    frequency1 = 0.5

    omega1 = 2 * np.pi * frequency1

    # Radial distance
    R = np.sqrt(X**2 + Y**2)

    # Radial inward direction (unit vectors)
    eps = 1e-8
    Ux_dir = -X / (R + eps)
    Uy_dir = -Y / (R + eps)

    # Two inward traveling radial waves
    t = 0
    phi1 = k1 * R + omega1 * t

    A = np.cos(phi1)

    # Final vector field
    U = A * Ux_dir
    V = A * Uy_dir

    U *= 0.0001
    V *= 0.0001

    return U, V


_globals = dict(ts=0)


def update(postOp):
    from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim

    live = postOp.live
    sim = live.simulation

    _globals["ts"] += 1
    ts = _globals["ts"]

    if ts % 100 != 0:
        return
    # print("ts", ts)

    hier = None
    hier = hierarchy_from_sim(live, qty="particles", pop="protons")
    L0 = hier.level(0, hier.times()[0])
    vmax = 0
    for ip, patch in enumerate(L0.patches):
        for i, name in enumerate(patch.patch_datas.keys()):
            pd = patch.patch_datas[name]
            # print("v", pd.dataset[0].v)
            for i in range(pd.dataset.size()):
                p = pd.dataset[i]
                # print("ptype", type(p))
                vmax = max(vmax, max(p.v))

    print("\nVmax: ", vmax)

    return
    if ts % 1000 != 0:
        return
    # print("ts++", ts)

    hier = None
    for i, c in enumerate(["x", "y"]):
        hier = hierarchy_from_sim(live, qty=f"EM_B_{c}", hier=hier)
    # for lvl_nbr, level in hier.levels(hier.times()[0]).items():
    L0 = hier.level(0, hier.times()[0])
    for ip, patch in enumerate(L0.patches):
        for i, name in enumerate(patch.patch_datas.keys()):
            pd = patch.patch_datas[name]
            nbrGhosts = pd.ghosts_nbr
            select = tuple([slice(nbrGhost, -(nbrGhost)) for nbrGhost in nbrGhosts])
            pd[pd.box] += b2(sim, *pd.meshgrid())[i][select]


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
            "timestamps": [final_time],
            "restart_time": "auto",
            "keep_last": 5,
        },
        strict="very",
    )

    def density(x, y):
        return 0.5

    def b(x, y):
        return b2(sim, x, y)

    def bx(x, y):
        return b(x, y)[0]

    def by(x, y):
        return b(x, y)[1]

    def bz(x, y):
        return 0

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return 0.00001

    C = "xyz"
    vvv = {
        **{f"vbulk{c}": vxyz for c in C},
        **{f"vth{c}": vthxyz for c in C},
        "nbr_part_per_cell": ppc,
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
