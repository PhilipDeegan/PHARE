#
#
#


import sys
import numpy as np

import pyphare.pharein as ph


cells = (512,)
dl = (1.0 / cells[0],)
L = 1.0
x0 = L / 2
sigma = 2.0 / 128

v_fast = np.sqrt(5.0 / 3.0 + 1.0)
time_step = 0.8 * dl[0] / v_fast
final_time = 1.0
timestamps = np.linspace(0, final_time, 41)


def density(x):
    return 1.0


def vx(x):
    return 0.08 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def vy(x):
    return 0.0


def vz(x):
    return 0.0


def bx(x):
    return 1.0


def by(x):
    return 0.0


def bz(x):
    return 0.0


def p(x):
    return 1.0


def config(*, diagdir):
    sim = ph.Simulation(
        cells=cells,
        dl=dl,
        time_step=time_step,
        final_time=final_time,
        boundary_types="periodic",
        refinement="tagging",
        max_nbr_levels=1,
        max_mhd_level=1,
        interp_order=2,
        gamma=5.0 / 3.0,
        mhd_timestepper="TVDRK3",
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        model_options=["MHDModel"],
        diag_options={
            "format": "phareh5",
            "options": {"dir": diagdir, "mode": "overwrite"},
        },
    )

    ph.MHDModel(
        density=density,
        vx=vx,
        vy=vy,
        vz=vz,
        bx=bx,
        by=by,
        bz=bz,
        p=p,
    )

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    return sim


def main():
    from pyphare.simulator.simulator import Simulator

    diagdir = sys.argv[1]
    Simulator(config(diagdir=diagdir), print_one_line=True).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    main()
