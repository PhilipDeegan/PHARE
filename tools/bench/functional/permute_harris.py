#
#
#


import os
import sys
import shutil
import datetime
import itertools
import subprocess
import numpy as np
from pathlib import Path
from copy import deepcopy

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.simulator.monitoring import MonitoringOptions
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI


def datetime_now():
    return datetime.datetime.now().replace(microsecond=0).strftime("%Y%m%d-%H%M%S")


MAKE_TAR_FILE = True
log_dir = Path(".log")
local_dir = Path(".phare")
out_dir = Path(f".phare_bench/{datetime_now()}")
monitoring_options = MonitoringOptions(interval=5, rank_modulo=1)

### test defaults
ndim = 3
dl = [0.4] * ndim
time_step = 0.001
final_time = 0.001
cells = [30, 30, 30]
ppc = 100


def _cells():
    ret = []
    for by in [1]:
        ret.append(deepcopy(cells))
        ret[-1][0] *= by
    return ret


permutables = [
    ("interp_order", [1]),
    ("cells", _cells()),
    ("ppc", [100]),
    ("max_nbr_levels", [2]),
    ("tag_buffer", [3, 4, 5]),
    ("tagging_threshold", [0.1, 0.2, 0.3, 0.4]),
]


def make_tarfile(outdir, outfile=None):
    if MAKE_TAR_FILE and cpp.mpi_rank() == 0:
        path = outdir.parent
        outfile = outfile or outdir.name
        args = ["tar", "czf", f"{path/outfile}.tar.gz", "-C", str(path), outdir.name]
        ec = subprocess.call(args)
        assert ec == 0, f"{ec}"
        clean_dir(outdir)


def mkdir(dir):
    dir.mkdir(parents=True, exist_ok=True)


def clean_dir(dir):
    if dir.exists():
        shutil.rmtree(str(dir))


def config(**kwargs):
    L = 0.5

    sim = ph.Simulation(
        refinement="tagging",
        interp_order=kwargs.get("interp_order", 1),
        time_step=kwargs.get("time_step", time_step),
        final_time=kwargs.get("final_time", final_time),
        dl=kwargs.get("dl", dl),
        cells=kwargs.get("cells", cells),
        max_nbr_levels=kwargs.get("max_nbr_levels", 1),
        tag_buffer=kwargs.get("tag_buffer", 3),
        nesting_buffer=kwargs.get("nesting_buffer", 1),
        tagging_threshold=kwargs.get("tagging_threshold", 0.1),
        hyper_resistivity=kwargs.get("hyper_resistivity", 0.008),
        hyper_mode=kwargs.get("hyper_mode", "spatial"),
        resistivity=kwargs.get("resistivity", 0.001),
        write_reports=False,
        strict=False,
    )

    def density(x, y, z):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y, z):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        return dBy1 + dBy2

    def bx(x, y, z):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y, z):
        return 0.0

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    def T(x, y, z):
        K = 0.7
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y, z):
        return 0.0

    def vy(x, y, z):
        return 0.0

    def vz(x, y, z):
        return 0.0

    def vthx(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthy(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthz(x, y, z):
        return np.sqrt(T(x, y, z))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": kwargs.get("ppc", ppc),
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "init": {"seed": cpp.mpi_rank() + 12},
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    return sim


def post_sim(summary, kwargs):
    log = local_dir / "summary.txt"
    with open(log, "w") as f:
        f.write(f"{kwargs}\n")
        f.write(summary)
    to = out_dir / datetime_now()
    shutil.copytree(log_dir, to / log_dir.name)
    shutil.copytree(local_dir, to / local_dir.name)
    make_tarfile(to)


def execute_permutation(keys, permutation):
    ph.global_vars.sim = None
    kwargs = {keys[i]: permutation[i] for i in range(len(keys))}
    sim = Simulator(config(**kwargs)).run(monitoring=monitoring_options)
    summary = sim.summary()
    sim.reset()
    cpp.mpi_barrier()
    if cpp.mpi_rank() == 0:
        post_sim(summary, kwargs)
    cpp.mpi_barrier()


### following function is called during test_case generation ###
def execute():
    if cpp.mpi_rank() == 0:
        mkdir(out_dir)  # should be unique
        mkdir(log_dir)

    cpp.mpi_barrier()
    keys = tuple(e[0] for e in permutables)
    for permutation in itertools.product(*(e[1] for e in permutables)):
        execute_permutation(keys, permutation)
    make_tarfile(out_dir)


if __name__ == "__main__":
    startMPI()
    execute()
