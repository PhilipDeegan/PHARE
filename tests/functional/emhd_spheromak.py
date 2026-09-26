#!/usr/bin/env python3
"""
Hall-only (EMHD) limit of the hybrid loop: ions frozen by a huge mass, cold,
uniform n, Te=0. Then Ve = -J/n, E = J x B / n, dB/dt = -curl(J x B / n):
E.J = 0 so int B^2 is conserved, and force-free fields (J || B) are exact
equilibria.

IC: Chandrasekhar-Kendall spheromak (curl B = k B, kR = 4.4934) inside r<R,
matched to a uniform anti-aligned field + dipole outside (J = 0), B continuous
at r=R so there is no surface current. PERTURB scales the toroidal component
(stays div-free, breaks J || B).

    python3 emhd_spheromak.py          # run + plot
    PLOT_ONLY=1 python3 emhd_spheromak.py
    PERTURB=0.2 python3 emhd_spheromak.py
"""
import os
import sys
import numpy as np
from pathlib import Path
from scipy.special import spherical_jn

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

PERTURB = float(os.environ.get("PERTURB", "0"))
ION_MASS = float(os.environ.get("ION_MASS", "1e6"))
NU = float(os.environ.get("NU", "1e-4"))
ETA = float(os.environ.get("ETA", "0"))

cells = (64, 64, 64)
dl = (0.125, 0.125, 0.125)
L = np.array(cells) * np.array(dl)
R = 1.5
B0 = 1.0
X0 = 4.493409457909064
K = X0 / R
J0_X0 = np.sin(X0) / X0
B_EXT = B0 * J0_X0 / 1.5

time_step = 2.5e-4
final_time = 10.0
n_dumps = 100
timestamps = np.linspace(0, final_time, n_dumps + 1)

tag = f"emhd_spheromak_p{PERTURB:g}"
diag_outputs = f"phare_outputs/test/{tag}"


def spheromak(x, y, z):
    dx, dy, dz = x - L[0] / 2, y - L[1] / 2, z - L[2] / 2
    rc2 = dx**2 + dy**2
    r = np.sqrt(rc2 + dz**2)
    rc = np.sqrt(rc2)
    rs = np.maximum(r, 1e-12)
    ct, st = dz / rs, rc / rs
    cp = np.where(rc > 0, dx / np.maximum(rc, 1e-12), 1.0)
    sp = np.where(rc > 0, dy / np.maximum(rc, 1e-12), 0.0)

    kr = K * rs
    j1 = spherical_jn(1, kr)
    j1p = spherical_jn(1, kr, derivative=True)
    j1_over = np.where(kr > 1e-6, j1 / np.maximum(kr, 1e-12), 1.0 / 3.0)
    Br_in = 2 * B0 * j1_over * ct
    Bt_in = -B0 * st * (j1_over + j1p)
    Bp_in = B0 * j1 * st * (1 + PERTURB)

    q = (R / rs) ** 3
    Br_out = B_EXT * ct * (1 - q)
    Bt_out = -B_EXT * st * (1 + q / 2)

    inside = r < R
    Br = np.where(inside, Br_in, Br_out)
    Bt = np.where(inside, Bt_in, Bt_out)
    Bp = np.where(inside, Bp_in, 0.0)

    Brho = Br * st + Bt * ct
    bz = Br * ct - Bt * st
    bx = Brho * cp - Bp * sp
    by = Brho * sp + Bp * cp
    return bx, by, bz


def config():
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        dl=dl,
        cells=cells,
        hyper_resistivity=NU,
        resistivity=ETA,
        diag_options={"format": "pharevtkhdf", "options": {"dir": diag_outputs}},
        strict=False,
    )

    C = "xyz"
    protons = {
        "charge": 1,
        "mass": ION_MASS,
        "density": lambda x, y, z: 1.0,
        **{f"vbulk{c}": (lambda x, y, z: 0.0) for c in C},
        **{f"vth{c}": (lambda x, y, z: 1e-6) for c in C},
        "nbr_part_per_cell": 50,
        "init": {"seed": 12334},
    }
    ph.MaxwellianFluidModel(
        bx=lambda x, y, z: spheromak(x, y, z)[0],
        by=lambda x, y, z: spheromak(x, y, z)[1],
        bz=lambda x, y, z: spheromak(x, y, z)[2],
        protons=protons,
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=timestamps)
    ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=timestamps)
    return sim


def global_array(vf, comp):
    shape = None
    out = None
    for patch in vf.level(0).patches:
        pd = patch.patch_datas[comp]
        data = np.asarray(pd[pd.box])
        if out is None:
            shape = tuple(np.array(cells) + 1)
            out = np.zeros(shape)
        lo = pd.box.lower
        out[tuple(slice(l, l + s) for l, s in zip(lo, data.shape))] = data
    return out


def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    run = Run(diag_outputs)
    plot_dir = Path(f"{diag_outputs}_plots")
    plot_dir.mkdir(parents=True, exist_ok=True)
    times = sorted(float(t) for t in run.all_times()["EM_B"])
    dV = float(np.prod(dl))
    mid = cells[1] // 2
    ext = [0, L[0], 0, L[2]]
    hist = dict(t=[], WB=[], fxb=[], vmax=[])
    for t in times:
        B = run.GetB(t)
        bx, by, bz = (global_array(B, c) for c in "xyz")
        hist["t"].append(t)
        hist["WB"].append(0.5 * float((bx**2 + by**2 + bz**2)[:-1, :-1, :-1].sum()) * dV)
        jx = np.gradient(bz, dl[1], axis=1) - np.gradient(by, dl[2], axis=2)
        jy = np.gradient(bx, dl[2], axis=2) - np.gradient(bz, dl[0], axis=0)
        jz = np.gradient(by, dl[0], axis=0) - np.gradient(bx, dl[1], axis=1)
        fx, fy, fz = jy * bz - jz * by, jz * bx - jx * bz, jx * by - jy * bx
        jb = np.sqrt(jx**2 + jy**2 + jz**2) * np.sqrt(bx**2 + by**2 + bz**2)
        hist["fxb"].append(float(np.sqrt(fx**2 + fy**2 + fz**2).sum() / max(jb.sum(), 1e-30)))
        try:
            V = run.GetVi(t)
            hist["vmax"].append(max(float(np.abs(global_array(V, c)).max()) for c in "xyz"))
        except Exception:
            hist["vmax"].append(np.nan)

        fig, ax = plt.subplots(2, 2, figsize=(12, 11))
        for a, v, name in (
            (ax[0, 0], by[:, mid, :], "B_y (toroidal at y=mid)"),
            (ax[0, 1], bz[:, mid, :], "B_z"),
            (ax[1, 0], jy[:, mid, :], "J_y"),
        ):
            m = max(float(np.abs(v).max()), 1e-30)
            im = a.imshow(v.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m)
            a.add_patch(plt.Circle((L[0] / 2, L[2] / 2), R, fill=False, ls=":", color="k"))
            fig.colorbar(im, ax=a, shrink=0.8)
            a.set_title(f"{name}, max |.|={m:.3e}")
            a.set_xlabel("x"); a.set_ylabel("z")
        a = ax[1, 1]
        th = np.array(hist["t"])
        a.plot(th, np.array(hist["WB"]) / hist["WB"][0], label="W_B / W_B(0)")
        a.plot(th, hist["fxb"], label="<|JxB|> / <|J||B|>")
        a.set_xlim(0, final_time); a.set_xlabel("t"); a.legend(fontsize=8)
        a2 = a.twinx()
        a2.semilogy(th, np.maximum(hist["vmax"], 1e-30), "k:", lw=0.8)
        a2.set_ylabel("max |V_i| (dotted)")
        fig.suptitle(f"{tag}  t={t:.3f}  ion mass={ION_MASS:g}  nu={NU:g} eta={ETA:g}")
        fig.tight_layout()
        fig.savefig(plot_dir / f"frame_{t:08.3f}.png", dpi=80)
        plt.close(fig)
    np.savez(plot_dir / "history.npz", **{k: np.array(v) for k, v in hist.items()})
    print(f"W_B drift {hist['WB'][-1] / hist['WB'][0] - 1:+.3e}, "
          f"force-free residual {hist['fxb'][0]:.3e} -> {hist['fxb'][-1]:.3e}")


if ph.PHARE_EXE:
    config()
elif __name__ == "__main__":
    startMPI()
    if not os.environ.get("PLOT_ONLY"):
        Simulator(config()).run()
    if cpp.mpi_rank() == 0:
        plot()
