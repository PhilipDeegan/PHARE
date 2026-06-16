#!/usr/bin/env python3
import os
import sys
import numpy as np
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

np.set_printoptions(threshold=sys.maxsize)

os.environ["PHARE_SCOPE_TIMING"] = "0"  # turn on scope timing


ph.NO_GUI()


cells = (101, 101, 101)
# cells = (41, 41, 41)
dl = (0.1, 0.1, 0.1)
dx, dy, dz = dl

name = "plasmoid"
diag_outputs = f"phare_outputs/test/{name}"

time_step = 0.00001
final_time = time_step  # 10
timestamps = [final_time]
# timestamps = np.arange(0, final_time + time_step, final_time / 10)
print("timestamps", timestamps)

params = dict(
    knot_radius=0.5,
    B0=5.0,
)


def b3(sim, x, y, z):
    # 1. Coordinate Centering (The Shift)
    # Extract domain and calculate center (midpoint)
    L = sim.simulation_domain()
    mid = np.array(L) / 2
    x0, y0, z0 = mid[0], mid[1], mid[2]

    # Local coordinate system relative to the center of the structure
    dx = x - x0
    dy = y - y0
    dz = z - z0

    # 2. Physical Scale of the "Bound Pressure System"
    # knot_radius: The geometric boundary of the proton
    # In PHARE tests, 'params' is often an attribute of sim or passed globally
    a = params.get("knot_radius", 1.0)
    # k: Taylor state constant for 3D stability (4.493 / a)
    k = 4.493 / a

    # 3. Distance Calculations in the Local Frame
    r_sq = dx**2 + dy**2 + dz**2
    r = np.sqrt(r_sq)
    R = np.sqrt(dx**2 + dy**2)  # Local cylindrical radius

    # Avoid division by zero at the center (the Bore)
    eps = 1e-10
    r_safe = np.where(r == 0, eps, r)
    R_safe = np.where(R == 0, eps, R)

    # 4. Spheromak Components (Internal Pressure & Twist)
    B0 = params.get("B0", 1.0)
    kr = k * r_safe
    sin_kr = np.sin(kr)
    cos_kr = np.cos(kr)

    # pol_factor: Defines the loops (The "Skeleton" of the Proton)
    pol_factor = B0 * (sin_kr / (kr * r_safe) - cos_kr / r_safe)
    # tor_factor: Defines the twist (The "Centripetal" Tension)
    tor_factor = B0 * sin_kr / r_safe

    # 5. Mapping back to Cartesian Vector Field
    # Bx and By combine the loops and the toroidal "twist"
    bx = pol_factor * (dx * dz / r_safe) - tor_factor * (dy / R_safe)
    by = pol_factor * (dy * dz / r_safe) + tor_factor * (dx / R_safe)

    # Bz creates the vertical core (The Bore)
    bz = pol_factor * (-(dx**2 + dy**2) / r_safe)

    # 6. Containment Mask
    # The plasmoid is a localized protrusion in the substrate
    mask = np.where(r < a, 1.0, 0.0)

    return bx * mask, by * mask, bz * mask


_globals = dict(ts=0)


def update(postOp):
    from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim

    live = postOp.live
    sim = live.simulation

    _globals["ts"] += 1
    ts = _globals["ts"]

    return
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
            # "dir": "checkpoints",
            "mode": "overwrite",
            # "elapsed_timestamps": [0],
            # "timestamps": [final_time],
            "restart_time": "auto",
        },
        strict=False,
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


if ph.PHARE_EXE:
    config()
elif __name__ == "__main__":
    Simulator(config()).run()
