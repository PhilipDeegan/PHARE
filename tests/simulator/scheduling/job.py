#!/usr/bin/env python3

import pyphare.pharein as ph

ph.Simulation(
    time_step_nbr=1,  # number of time steps (not specified if time_step and final_time provided)
    time_step=0.010,  # simulation final time (not specified if time_step and time_step_nbr provided)
    boundary_types=["periodic", "periodic"],
    cells=[100, 100],  # integer or tuple length == dimension
    dl=[0.1, 0.1],  # mesh size of the root level, float or tuple
)
density = lambda x: 2.0
bx, by, bz = (lambda x: 1 for i in range(3))
ex, ey, ez = (lambda x: 1 for i in range(3))
vx, vy, vz = (lambda x: 1.0 for i in range(3))
vthx, vthy, vthz = (lambda x: 1.0 for i in range(3))
vvv = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthx,
    "vthy": vthy,
    "vthz": vthz,
}

ph.MaxwellianFluidModel(
    bx=bx,
    by=by,
    bz=bz,
    protons={"charge": 1, "density": density, **vvv, "init": {"seed": 1337}},
)
ph.ElectronModel(closure="isothermal", Te=0.12)
