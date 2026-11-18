"""Tagger (pyphare.simulator.tagger) against a run: on a 2D Harris sheet with two levels, the
level 1 a run makes at initialization and the one Tagger finds from the initial fields, outside
the simulation, cover the same cells, for each tagging method."""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.data.wrangler import DataWrangler
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.tagger import Tagger
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = (200, 100)
dl = (0.40, 0.40)
min_overlap = 0.95
thresholds = {
    "default": 0.1,
    "lohner-per-direction": 0.3,
    "lohner-cross-derivatives": 0.3,
    "gradient-scale": {"threshold": 0.1, "ref_magnitude": 0.1},
}


def harris(tagging):
    sim = ph.Simulation(
        time_step=0.005,
        time_step_nbr=1,
        cells=cells,
        dl=dl,
        interp_order=1,
        refinement="tagging",
        max_nbr_levels=2,
        tagging=tagging,
        hyper_resistivity=0.002,
        resistivity=0.001,
    )
    Lx, Ly = sim.simulation_domain()
    L = 0.5

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def density(x, y):
        return (
            0.2
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def perturbation(x, y, y0):
        return 0.1 * np.exp(-((x - 0.5 * Lx) ** 2 + (y - y0) ** 2))

    def bx(x, y):
        y1, y2 = y - 0.3 * Ly, y - 0.7 * Ly
        dbx = -2 * y1 * perturbation(x, y, 0.3 * Ly) + 2 * y2 * perturbation(
            x, y, 0.7 * Ly
        )
        return -1 + 2 * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dbx

    def by(x, y):
        x0 = x - 0.5 * Lx
        return 2 * x0 * (perturbation(x, y, 0.3 * Ly) - perturbation(x, y, 0.7 * Ly))

    def bz(x, y):
        return 0.0 * x

    def T(x, y):
        b2 = bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2
        return np.maximum(1.0 / density(x, y) * (1 + 0.7 - 0.5 * b2), 0.1)

    def vth(x, y):
        return np.sqrt(T(x, y))

    def zero(x, y):
        return 0.0 * x

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 10,
            "init": {"seed": 12334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    return sim


def cells_of(boxes):
    """the set of cells of the boxes, as tuples"""
    out = set()
    for lower, upper in boxes:
        lower, upper = np.asarray(lower), np.asarray(upper)
        out.update(tuple(lower + i) for i in np.ndindex(*(upper - lower + 1)))
    return out


class TaggerTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None
        super().tearDown()

    def test_level_one_overlaps_the_run(self):
        for method, threshold in thresholds.items():
            with self.subTest(method=method):
                ph.global_vars.sim = None
                sim = harris({"method": method, "quantities": {"B": threshold}})

                tagged = Tagger(sim).tag()
                self.assertIn(1, tagged, "the Tagger refines nothing")
                tag_cells = cells_of((b.lower, b.upper) for b in tagged[1])

                self.simulator = Simulator(sim)
                self.simulator.initialize()
                wrangler = DataWrangler(self.simulator)
                self.assertEqual(wrangler.getNumberOfLevels(), 2)
                run_cells = cells_of(
                    (np.asarray(pd.lower), np.asarray(pd.upper))
                    for pd in wrangler.getPatchLevel(1).getB("x")
                )
                wrangler.kill()
                self.simulator.reset()
                self.simulator = None

                overlap = len(tag_cells & run_cells) / len(tag_cells | run_cells)
                print(f"{method}: level 1 overlap {overlap:.4f}")
                self.assertGreaterEqual(overlap, min_overlap)


if __name__ == "__main__":
    unittest.main()
