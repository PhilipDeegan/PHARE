#
#
#


import sys
import numpy as np

import pyphare.pharein as ph


def config(*, diagdir: str, T: float, final_time: float):
    cells = (512,)
    dl = (1.0 / cells[0],)
    L = 1.0

    v_fast = np.sqrt(5.0 / 3.0 * T + 1.0)
    time_step = 0.8 * dl[0] / v_fast

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
        density=lambda x: 1.0,
        vx=lambda x: np.sin(2 * np.pi / L * x) * 1.5,
        vy=lambda x: 0.0,
        vz=lambda x: 0.0,
        bx=lambda x: 1.0,
        by=lambda x: 0.0,
        bz=lambda x: 0.0,
        p=lambda x: T,
    )

    timestamps = np.linspace(0, final_time, 17)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    return sim


def main():
    from pyphare.simulator.simulator import Simulator

    if len(sys.argv) != 4:
        print('This code needs 3 parameters: diagdir, T, final_time')
        sys.exit(1)

    diagdir = sys.argv[1]
    T = float(sys.argv[2])
    final_time = float(sys.argv[3])

    Simulator(config(diagdir=diagdir, T=T, final_time=final_time), print_one_line=True).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    main()
