#!/usr/bin/env python3


import unittest
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

import numpy as np
import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib

cpp = cpp_lib()
startMPI()

mpi_size = cpp.mpi_size()
time_step_nbr = 5
time_step = 0.001
cells = (100, 100)
dl = (0.2, 0.2)
diag_outputs = "phare_outputs/harris/2d/load_balancing"
timestamps = [x * time_step for x in range(time_step_nbr + 1)]


def config():
    sim = ph.Simulation(
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        advanced={
            "integrator/rebalance_coarsest": 1,
            "integrator/rebalance_coarsest_every": 1,
            "integrator/rebalance_coarsest_on_init": 1,
            "integrator/flexible_load_tolerance": 0.05,
        },
        loadbalancing="nppc",
    )

    def ppc_by_icell(x, y):
        # print("ppc_by_icell(icell)", x)
        # print("ppc_by_icell(icell)", x.shape, y.shape)
        ppc = y.copy()
        ppc[:] = 90
        # print("ppc", ppc)
        ppc[np.where(np.isclose(y, 10, atol=0.1))] = 590
        # print("ppc", ppc)
        return 100  # ppc

    def density(x, y):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 1
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vxyz,
        "vbulky": vxyz,
        "vbulkz": vxyz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
        "nbr_part_per_cell": (ppc_by_icell, 100),
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.ParticleDiagnostics(
        quantity="domain", write_timestamps=timestamps, population_name="protons"
    )
    return sim


def get_time(path, time=0, datahier=None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    return hierarchy_from(
        h5_filename=path + "/ions_pop_protons_domain.h5", time=time, hier=datahier
    )


def time_info(sim, time=0):
    hier = get_time(diag_outputs, time)
    per_rank = {f"p{rank}": 0 for rank in range(mpi_size)}

    def _parse_rank(patch_id):
        return patch_id.split("#")[0]

    for ilvl, lvl in hier.levels().items():
        for patch in lvl:
            for pd_key, pd in patch.patch_datas.items():
                per_rank[_parse_rank(patch.id)] += pd.size()

    return per_rank


class LoadBalancingTest(SimulatorTest):
    def test_has_balanced(self):
        if mpi_size == 1:  # doesn't make sense
            return

        sim = config()
        self.register_diag_dir_for_cleanup(diag_outputs)
        Simulator(sim).run()

        if cpp.mpi_rank() > 0:
            return

        t0_sdev = np.std(list(time_info(sim).values()))
        tend_sdev = np.std(list(time_info(sim, timestamps[-1]).values()))

        assert tend_sdev < t0_sdev * 0.1  # empirical


if __name__ == "__main__":
    unittest.main()
