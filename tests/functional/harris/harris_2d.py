#!/usr/bin/env python3



import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI


from tests.simulator import SimulatorTest


ph.NO_GUI()



cpp = cpp_lib()
diag_outputs = "phare_outputs/test/harris/2d"
time_step_nbr = 1000
time_step = 0.001
final_time = time_step * time_step_nbr
cells = (200, 400)


def default_timestamps():
    dt = 10 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)
    return timestamps


def default_setup():
    return ph.Simulation(
        # smallest_patch_size=15,
        # largest_patch_size=25,
        # time_step_nbr=time_step_nbr,
>>>>>>> 9a9e694f (realtime plots)
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
    )


def config(sim=None, timestamps=None, seed=12334):
    L = 0.5
    if sim is None:
        sim = default_setup()
    if timestamps is None:
        timestamps = default_timestamps()

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
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

    def bx(x, y):
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

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
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

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": seed}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["mass_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["density", "pressure_tensor"]:
        ph.FluidDiagnostics(
            quantity=quantity, write_timestamps=timestamps, population_name="protons"
        )

    ph.InfoDiagnostics(quantity="particle_count")
    ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

    return sim



def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    pop_name = "protons"
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        )
        run.GetN(time, pop_name=pop_name).plot(
            filename=plot_file_for_qty(plot_dir, "N", time), plot_patches=True
        )
        for c in ["x", "y", "z"]:
            run.GetB(time).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetJ(time).plot(
            filename=plot_file_for_qty(plot_dir, "jz", time),
            qty="z",
            plot_patches=True,
            vmin=-2,
            vmax=2,
        )
        run.GetPressure(time, pop_name=pop_name).plot(
            filename=plot_file_for_qty(plot_dir, "Pxx", time),
            qty=pop_name + "_Pxx",
            plot_patches=True,
            vmin=0,
            vmax=2.7,
        )
        run.GetPressure(time, pop_name=pop_name).plot(
            filename=plot_file_for_qty(plot_dir, "Pzz", time),
            qty=pop_name + "_Pzz",
            plot_patches=True,
            vmin=0,
            vmax=1.5,
        )


def realtime_plots(new_time, post_op):
    if cpp.mpi_rank() == 0:
        from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim
        import matplotlib.pyplot as plt

        plt.clf()

        axes = []

        def add_ax():
            axes.append(
                post_op.fig.add_subplot(post_op.rows, post_op.cols, len(axes) + 1)
            )
            return axes[-1]

        for c in ["x", "y", "z"]:
            hierarchy_from_sim(post_op.live, qty=f"EM_B_{c}").plot(
                ax=add_ax(), qty=f"{c}", plot_patches=True
            )

        plt.draw()
        plt.pause(0.01)


def b3(sim, x, y, z, lx, ly, lz):
    L = sim.simulation_domain()[0]
    mid = L / 2

    X = (x - mid).reshape(lx, ly, lz)
    Y = (y - mid).reshape(lx, ly, lz)
    Z = (z - mid).reshape(lx, ly, lz)

    U, V, W = -X, -Y, -Z

    # Normalize vectors to unit length
    magnitude = np.sqrt(U**2 + V**2 + W**2) + 1e-5  # Avoid division by zero

    # Normalize vectors (unit length)
    U /= magnitude
    V /= magnitude
    W /= magnitude

    # Define circular mask (radius = 25)
    radius = L * 0.4
    diff = 0.2

    # outer mask
    mask = X**2 + Y**2 + Z**2 <= (radius + diff) ** 2
    U[~mask] = 0
    V[~mask] = 0
    W[~mask] = 0

    # inner mask
    mask = X**2 + Y**2 + Z**2 >= (radius - diff) ** 2
    U[~mask] = 0
    V[~mask] = 0
    W[~mask] = 0

    U *= 0.001
    V *= 0.001
    W *= 0.001

    return U, V, W


def b2(sim, x, y, lx, ly):
    L = sim.simulation_domain()[0]
    mid = L / 2

    X = x - mid
    Y = y - mid

    U, V = -X, -Y

    # Normalize vectors to unit length
    magnitude = np.sqrt(U**2 + V**2) + 1e-5  # Avoid division by zero

    # Normalize vectors (unit length)
    U /= magnitude
    V /= magnitude

    # Define circular mask (radius = 25)
    radius = L * 0.4
    diff = 0.2

    # outer mask
    mask = X**2 + Y**2 <= (radius + diff) ** 2
    U[~mask] = 0
    V[~mask] = 0

    # inner mask
    mask = X**2 + Y**2 >= (radius - diff) ** 2
    U[~mask] = 0
    V[~mask] = 0

    U *= 0.001
    V *= 0.001

    return U, V


def update(post_op):
    from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim

    live = post_op.live
    sim = live.simulation

    hier = None
    for c in ["x", "y"]:
        hier = hierarchy_from_sim(live, qty=f"EM_B_{c}", hier=hier)
    for lvl_nbr, level in hier.levels(hier.times()[0]).items():
        for ip, patch in enumerate(level.patches):
            pdata_names = list(patch.patch_datas.keys())
            for i, name in enumerate(pdata_names):
                pd = patch.patch_datas[name]
                nbrGhosts = pd.ghosts_nbr
                select = tuple([slice(nbrGhost, -(nbrGhost)) for nbrGhost in nbrGhosts])
                pd[pd.box] += b2(sim, *pd.meshgrid(), *pd.size)[i][select]


class PostOp:
    def __init__(self, live, diag_dir):
        import matplotlib.pyplot as plt

        plt.ion()
        plt.show()

        self.live = live
        self.diag_dir = diag_dir
        self.cols = 3
        self.rows = 2
        self.scale = 4
        figsize = (self.cols * self.scale, self.rows * self.scale)
        self.fig = plt.figure(figsize=figsize)
        self.fig.tight_layout()

    def __call__(self, new_time):
        ...
        # realtime_plots(new_time, self)
        # update(self)




class HarrisTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(HarrisTest, self).__init__(*args, **kwargs)
        self.simulator = None


    def tearDown(self):
        super(HarrisTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        if cpp.mpi_rank() == 0:
            plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
            plot_dir.mkdir(parents=True, exist_ok=True)
            plot(diag_dir, plot_dir)
        cpp.mpi_barrier()
        return self


    # def test_run(self, diag_dir=None, sim=None):
    #     diag_dir = diag_dir if diag_dir else diag_outputs
    #     self.plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
    #     self.plot_dir.mkdir(parents=True, exist_ok=True)
    #     sim = sim if sim else config()
    #     self.register_diag_dir_for_cleanup(diag_dir)
    #     live = Simulator(sim)
    #     live.post_advance = PostOp(live, diag_dir)
    #     live.run().reset()
    #     return self

    # def plot(self, timestamps, diag_dir, plot_dir):
    #     run = self.getRun(diag_dir)
    #     for time in timestamps:
    #         run.GetDivB(time).plot(
    #             filename=plot_file_for_qty(plot_dir, "divb", time),
    #             plot_patches=True,
    #             vmin=1e-11,
    #             vmax=2e-10,
    #         )
    #         run.GetRanks(time).plot(
    #             filename=plot_file_for_qty(plot_dir, "Ranks", time),
    #             plot_patches=True,
    #         )
    #         run.GetN(time, pop_name="protons").plot(
    #             filename=plot_file_for_qty(plot_dir, "N", time),
    #             plot_patches=True,
    #         )
    #         for c in ["x", "y", "z"]:
    #             run.GetB(time).plot(
    #                 filename=plot_file_for_qty(plot_dir, f"b{c}", time),
    #                 qty=f"{c}",
    #                 plot_patches=True,
    #             )
    #         run.GetJ(time).plot(
    #             filename=plot_file_for_qty(plot_dir, "jz", time),
    #             qty="z",
    #             plot_patches=True,
    #             vmin=-2,
    #             vmax=2,
    #         )

    # def scope_timing(self, diag_dir):
    #     try:
    #         from tools.python3 import plotting as m_plotting

    #         m_plotting.plot_run_timer_data(diag_dir, cpp.mpi_rank())
    #     except ImportError:
    #         print("Phlop not found - install with: `pip install phlop`")
    #     except FileNotFoundError:
    #         print("Phlop installed but not active")


if __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()
