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

name = "plasmoids"
diag_outputs = f"phare_outputs/test/{name}"

time_step = 0.00001
final_time = 0.0001
restart_options = {
    "dir": "checkpoints",
    "mode": "overwrite",
    "timestamps": [final_time],
    "restart_time": "auto",
}
max_nbr_levels = 1
start_time = ph.restarts.restart_time(restart_options) or 0
timestamps = [start_time, final_time]
# timestamps = np.arange(0, final_time + time_step, final_time / 10)
print("timestamps", timestamps)

params = dict(
    knot_radius=0.5,
    B0=5.0,
)


def get_spheromak(x, y, z, cx, cy, cz, rotation=1.0):
    # Local Shift
    dx, dy, dz = x - cx, y - cy, z - cz
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    R = np.sqrt(dx**2 + dy**2)

    eps = 1e-10
    r_safe = np.where(r == 0, eps, r)
    R_safe = np.where(R == 0, eps, R)

    a = params.get("knot_radius", 1.0)
    k = 4.493 / a
    B0 = params.get("B0", 1.0)

    kr = k * r_safe
    pol_f = B0 * (np.sin(kr) / (kr * r_safe) - np.cos(kr) / r_safe)
    # 'rotation' allows us to test Same-Spin vs Opposite-Spin
    tor_f = B0 * np.sin(kr) / r_safe * rotation

    bx = pol_f * (dx * dz / r_safe) - tor_f * (dy / R_safe)
    by = pol_f * (dy * dz / r_safe) + tor_f * (dx / R_safe)
    bz = pol_f * (-(dx**2 + dy**2) / r_safe)

    mask = np.where(r < a, 1.0, 0.0)
    return bx * mask, by * mask, bz * mask


def b3(sim, x, y, z):
    L = sim.simulation_domain()
    mid = np.array(L) / 2
    d = params.get("separation", 1.0)

    # 1. Start with a Turbulent "Sea" (The Substrate)
    # Using a high-frequency hash or sine to simulate "Zero Point" pressure
    sea_intensity = 15
    bx_sea = sea_intensity * np.sin(10 * x) * np.cos(10 * y)
    by_sea = sea_intensity * np.sin(10 * y) * np.cos(10 * z)
    bz_sea = sea_intensity * np.sin(10 * z) * np.cos(10 * x)

    # 2. Layer the 5x5x5 "Loophole" Grid
    bx_knot, by_knot, bz_knot = 0, 0, 0

    # Iterate through a 5x5x5 grid centered on 'mid'
    # Range is -2 to 2 to give us 5 units (-2, -1, 0, 1, 2)
    lo, up = -4, 5
    for i in range(lo, up):
        for j in range(lo, up):
            for k in range(lo, up):
                cx = mid[0] + i * d
                cy = mid[1] + j * d
                cz = mid[2] + k * d

                # Superimpose each knot
                rotation = 1.0 if (i + j + k) % 2 == 0 else -1.0
                bx, by, bz = get_spheromak(x, y, z, cx, cy, cz, rotation=rotation)
                bx_knot += bx
                by_knot += by
                bz_knot += bz

    # 3. Combine: The knots "Organize" the sea
    # Where the knot exists, it replaces the random noise with coherent structure
    return bx_sea + bx_knot, by_sea + by_knot, bz_sea + bz_knot


def b3_125_a(sim, x, y, z):  # 125 grid
    L = sim.simulation_domain()
    mid = np.array(L) / 2

    d = params.get("separation", 2.0)

    bx_total, by_total, bz_total = 0, 0, 0

    # Iterate through a 5x5x5 grid centered on 'mid'
    # Range is -2 to 2 to give us 5 units (-2, -1, 0, 1, 2)
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                cx = mid[0] + i * d
                cy = mid[1] + j * d
                cz = mid[2] + k * d

                # Superimpose each knot
                rotation = 1.0 if (i + j + k) % 2 == 0 else -1.0
                bx, by, bz = get_spheromak(x, y, z, cx, cy, cz, rotation=rotation)
                bx_total += bx
                by_total += by
                bz_total += bz

    return bx_total, by_total, bz_total


def b3_6(sim, x, y, z):
    L = sim.simulation_domain()
    mid = np.array(L) / 2

    d = params.get("separation", 2.0)

    # Define the 6 directions for the 3D "+" sign
    offsets = [[d, 0, 0], [-d, 0, 0], [0, d, 0], [0, -d, 0], [0, 0, d], [0, 0, -d]]

    bx_total, by_total, bz_total = 0, 0, 0

    for off in offsets:
        cx, cy, cz = mid + np.array(off)
        # Call the single spheromak function for each center
        bx, by, bz = get_spheromak(x, y, z, cx, cy, cz)
        bx_total += bx
        by_total += by
        bz_total += bz

    return bx_total, by_total, bz_total


def b3_2(sim, x, y, z):
    L = sim.simulation_domain()
    mid = np.array(L) / 2

    # 1. Define the Offset for the two "Protons"
    # We'll place them on either side of the midpoint along the X-axis
    separation = params.get("separation", 2.0)

    # Centers for Knot A and Knot B
    center_a = mid.copy()
    center_a[0] -= separation / 2

    center_b = mid.copy()
    center_b[0] += separation / 2

    # 2. Generate both knots
    bxa, bya, bza = get_spheromak(*center_a, rotation=1.0)
    bxb, byb, bzb = get_spheromak(*center_b, rotation=1.0)

    # 3. Superimpose (Linear Sum)
    return bxa + bxb, bya + byb, bza + bzb


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
        max_nbr_levels=max_nbr_levels,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "pharevtkhdf",
            "options": {"dir": diag_outputs},
        },
        restart_options=restart_options,
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
