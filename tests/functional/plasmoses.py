#!/usr/bin/env python3
import os
import sys
import numpy as np
from pathlib import Path
from scipy.special import j0, j1, spherical_jn

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

np.set_printoptions(threshold=sys.maxsize)

os.environ["PHARE_SCOPE_TIMING"] = "0"  # turn on scope timing


ph.NO_GUI()


cells = (101, 101, 101)
mid = (51, 51, 51)
# cells = (41, 41, 41)
dl = (0.1, 0.1, 0.1)
dx, dy, dz = dl

diag_outputs = f"phare_outputs/test/plasmoses"

time_step = 0.001
final_time = 10
restart_options = {
    "dir": "checkpoints",
    "mode": "overwrite",
    "timestamps": [final_time],
    "restart_time": "auto",
}
max_nbr_levels = 1
start_time = ph.restarts.restart_time(restart_options) or 0
timestamps = [start_time, final_time]
timestamps = np.arange(0, final_time + time_step, final_time / 100)
print("timestamps", timestamps)

params = dict(
    knot_radius=0.5,
    B0=5.0,
)


def get_taylor_spheromak(x, y, z, cx, cy, cz, R=1.0, B0=1.0):
    # 1. Calculate relative coordinates and radial distance
    dx, dy, dz = x - cx, y - cy, z - cz
    r = np.sqrt(dx**2 + dy**2 + dz**2)

    # 2. Create the Mask (The "Neutral Zone" boundary)
    # This prevents the ValueError by handling the whole array at once
    mask = (r <= R) & (r > 0)

    # Initialize output arrays with zeros
    bx = np.zeros_like(x)
    by = np.zeros_like(y)
    bz = np.zeros_like(z)

    # 3. Calculate only for points inside the "Ball of Potential"
    k = 4.493 / R
    rho = k * r[mask]

    # Spherical Bessel j1(rho)
    j1_rho = spherical_jn(1, rho)

    # Force-Free projection factor
    # In a discrete substrate, this is the "Integer Gearing"
    common_factor = B0 * j1_rho / r[mask]

    # 4. Populate the vector components
    bx[mask] = common_factor * (dz[mask] * dx[mask] / r[mask] - dy[mask])
    by[mask] = common_factor * (dz[mask] * dy[mask] / r[mask] + dx[mask])
    bz[mask] = common_factor * (-(dx[mask] ** 2 + dy[mask] ** 2) / r[mask])

    return bx, by, bz


def get_taylor_spheromak_cylinder(x, y, z, cx, cy, cz, R=1.0, B0=1.0):
    # 1. Calculate local coordinates relative to center
    dx, dy, dz = x - cx, y - cy, z - cz
    r_cyl = np.sqrt(dx**2 + dy**2)
    r_sph = np.sqrt(dx**2 + dy**2 + dz**2)

    # 2. Define the boundary (The Neutral Zone)
    # The first zero of J1 is approx 3.8317
    k = 3.8317 / R

    if all(r_sph > R):
        return 0.0, 0.0, 0.0

    # 3. Woltjer-Taylor Force-Free Solution
    # This geometry creates a 'locked' ball of potential
    # B_z is poloidal, B_phi is toroidal

    arg = k * r_cyl
    # Poloidal component (vertical through the bore)
    bz = B0 * j0(arg)

    # Toroidal component (the 'twist' around the bore)
    # We use the cross product logic to get bx and by from b_phi
    b_phi = B0 * j1(arg)

    if all(r_cyl) == 0:
        bx, by = 0.0, 0.0
    else:
        bx = -b_phi * (dy / r_cyl)
        by = b_phi * (dx / r_cyl)

    return bx, by, bz


def b3(sim, x, y, z):
    L = sim.simulation_domain()
    mid = np.array(L) / 2

    cx = mid[0]
    cy = mid[1]
    cz = mid[2]

    bx, by, bz = get_taylor_spheromak(x, y, z, cx, cy, cz)

    return bx, by, bz


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

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            component_properties={
                "x": [
                    ph.FieldDiagnosticSlice(
                        (mid[0], 0, 0), (mid[0], cells[1], cells[2])
                    )
                ],
                # "y": [ph.FieldDiagnosticSlice((20, 0, 0), (20, 40, 40))],
                # "z": [ph.FieldDiagnosticSlice((20, 0, 0), (20, 40, 40))],
            },
        )

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/{qty}_t{time}.png"


def plot():
    import shlex
    import subprocess

    run = Run(diag_outputs)
    plot_dir = Path(f"{diag_outputs}_plots")
    plot_dir.mkdir(parents=True, exist_ok=True)

    times = run.all_times()["EM_B"]
    for time in times:
        B = run.GetB(time, all_primal=False)
        B = B.magnitude()
        B.finest().plot(filename=plot_file_for_qty(plot_dir, f"Bmag", time))

    cmd = shlex.split(
        f"ffmpeg -r 10 -y -pattern_type glob -i '{plot_dir}/Bmag_t*.png' -c:v libx264 -crf 0 Bmagg.mp4"
    )

    subprocess.call(cmd)


if ph.PHARE_EXE:
    config()
elif __name__ == "__main__":
    startMPI()
    # Simulator(config()).run()
    if cpp.mpi_rank() == 0:
        plot()
