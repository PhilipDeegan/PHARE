#!/usr/bin/env python3
import os
import numpy as np

from pathlib import Path
import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI
from tests.diagnostic import dump_all_diags

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing
"""
  For scope timings to work
  The env var PHARE_SCOPE_TIMING must be == "1" (or "true")
    See src/phare/phare.hpp
  CMake must be configured with: -DwithPhlop=ON
  And a LOG_LEVEL must be defined via compile args: -DPHARE_LOG_LEVEL=1
  Or change the default value in src/core/logger.hpp
  And phlop must be available on PYTHONPATH either from subprojects
   or install phlop via pip
"""


ph.NO_GUI()
cpp = cpp_lib()
startMPI()


diag_outputs = "phare_outputs/test/harris/2d"
time_step_nbr = 1
time_step = 0.001
final_time = time_step * time_step_nbr
dt = 10 * time_step
nt = final_time / dt + 1
timestamps = [0, final_time]  # dt * np.arange(nt)
cells = (12, 12)


def config(diag_dir):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        # smallest_patch_size=4,
        # largest_patch_size=25,
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=(0.2, 0.2),
        refinement="tagging",
        max_nbr_levels=1,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=False,
        tag_buffer=4,
    )

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

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 100,
    }

    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    dump_all_diags(model.populations)
    # for quantity in ["E", "B"]:
    #     ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    # ph.InfoDiagnostics(quantity="particle_count")  # defaults all coarse time steps

    # for quantity in ["density", "bulkVelocity"]:
    #     ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    # ph.FluidDiagnostics(
    #     quantity="density", write_timestamps=timestamps, population_name="protons"
    # )

    return sim


def get_time(path, time, datahier=None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", times=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", times=time, hier=datahier)
    return datahier


def cmp():
    ch = get_time(diag_outputs, 0)
    gh = get_time(diag_outputs + "_gpu", 0)
    assert ch == gh

    # for key in ch.level(0)[0].keys():
    #     print(ch.level(0)[0].box, gh.level(0)[0].box)
    #     assert ch.level(0)[0][key] == gh.level(0)[0][key], f"{key} failed"


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    Path(plot_dir).mkdir(exist_ok=True, parents=True)
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time),
            plot_patches=True,
        )
        run.GetN(time, pop_name="protons").plot(
            filename=plot_file_for_qty(plot_dir, "N", time),
            plot_patches=True,
        )
        for c in ["x", "y", "z"]:
            run.GetB(time).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                qty=f"{c}",
                plot_patches=True,
            )
        run.GetJ(time).plot(
            filename=plot_file_for_qty(plot_dir, "jz", time),
            qty="z",
            plot_patches=True,
            vmin=-2,
            vmax=2,
        )


def main():
    Simulator(config(diag_outputs)).setup().initialize().run()
    Simulator(config(diag_outputs + "_gpu")).setup(layout=3).initialize().run()

    try:
        from tools.python3 import plotting as m_plotting

        m_plotting.plot_run_timer_data(diag_outputs, cpp.mpi_rank())
    except ImportError:
        print("Phlop not found - install with: `pip install phlop`")

    except FileNotFoundError:
        print("Phlop installed but not active`")

    plot(diag_outputs, diag_outputs + "_plots")
    plot(diag_outputs + "_gpu", diag_outputs + "_gpu_plots")

    cmp()


if __name__ == "__main__":
    main()
