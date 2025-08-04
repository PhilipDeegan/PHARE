#
#
#

import os
import yaml
import inspect
import importlib
import numpy as np
from pathlib import Path
from pyphare.cpp import cpp_lib
from dataclasses import asdict as dc_asdict, dataclass, field as dc_field

cpp = cpp_lib()

USAGE = """!PHARE CLI!


example: python3 pyphare/pyphare --watch 1 tests.functional.harris.harris_2d.config

Executes Simulation object defined in config() function in file tests/functional/harris/harris_2d.py
            with realtime plots showing live data

"""


@dataclass
class CliHelp:
    config = "Enable cmake build config tests extraction"
    watch = "Interval to watch realtime plots, default 0 is disabled"
    monitor = "Interval to record system utility information per rank"


def cli_args_parser():
    import argparse

    _help = CliHelp()
    parser = argparse.ArgumentParser(
        description=USAGE, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config", default="", help=_help.config)
    parser.add_argument("--watch", default=0, help=_help.watch)
    parser.add_argument("--monitor", default=0, help=_help.monitor)
    parser.add_argument("remaining", nargs=argparse.REMAINDER)
    return parser


def verify_cli_args(cli_args):
    try:
        cli_args.watch = int(cli_args.watch)
        cli_args.monitor = int(cli_args.monitor)
    except ValueError:
        raise ValueError("Interval must be an integer")
    return cli_args


FORCE_RAISE_ON_IMPORT_ERROR = os.getenv(
    "PHLOP_FORCE_RAISE_ON_IMPORT_ERROR", "False"
).lower() in ("true", "1", "t")


def is_module(module, rethrow=False):
    try:
        importlib.import_module(module)
    except (ValueError, ModuleNotFoundError) as e:
        if rethrow:
            raise e
        return False
    return True


bad_input_error = (
    "bad input, must be module or module.function to execute to setup simulation"
)


def get_callable(remaining):
    if is_module(remaining):
        return remaining
    if not remaining.find("."):
        raise ValueError(bad_input_error)
    mod_bits = cli_args.remaining[0].split(".")
    path, fn = ".".join(mod_bits[:-1]), mod_bits[-1]
    if not is_module(path, rethrow=True):
        raise ValueError(bad_input_error)
    for name, el in inspect.getmembers(importlib.import_module(path)):
        if name == fn:
            return el
    raise ValueError(bad_input_error)


def realtime_plots(new_time, post_op):

    _locals = dict(idx=0)

    def add_or_get(i):
        if len(post_op.axes) <= i:
            post_op.axes.append(
                post_op.fig.add_subplot(post_op.rows, post_op.cols, i + 1)
            )
        return post_op.axes[i]

    def _try_plot_axis(qty, k):
        try:
            ax = add_or_get(_locals["idx"])
            hierarchy_from_sim(post_op.live, qty=qty).plot(
                ax=ax, qty=k, plot_patches=True
            )
            _locals["idx"] += 1
        except FileNotFoundError:
            ...

    def _try_plot_vecfield(qty):
        for c in ["x", "y", "z"]:
            _try_plot_axis(f"{qty}_{c}", c)

    if cpp.mpi_rank() == 0:
        from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim
        import matplotlib.pyplot as plt

        # plt.cla()

        _try_plot_vecfield("EM_B")
        _try_plot_vecfield("EM_E")
        _try_plot_vecfield("bulkVelocity")

        post_op.fig.canvas.draw_idle()
        plt.pause(0.05)




class PostOp:
    def __init__(self, cli_args, live):
        import matplotlib.pyplot as plt

        plt.ion()
        plt.show()
        self.live = live
        self.cols = 3
        self.rows = 3
        self.scale = 3
        figsize = (self.cols * self.scale, self.rows * self.scale)
        self.fig = plt.figure(figsize=figsize)
        self.fig.tight_layout()
        self.tidx = 0
        self.axes = []
        if cli_args.watch:
            realtime_plots(live.currentTime(), self)

    def __call__(self, new_time):
        if cli_args.watch and self.tidx % cli_args.watch == 0:
            realtime_plots(new_time, self)
        self.tidx += 1


def print_sim_summary(sim):
    if cpp.mpi_rank() > 0:
        return

    def _pr(s, tabs=0):
        print(f"{"\t" * tabs} {s}")

    print("Simulation summary:")
    _pr(f"cells: {sim.cells}", 1)
    _pr(f"final_time: {sim.final_time}", 1)
    _pr(f"time_step: {sim.time_step}", 1)
    _pr(f"time_step_nbr: {sim.time_step_nbr}", 1)
    _pr(f"lower bound memory requirement: todo", 1)


def read_yaml(cli_args):
    return yaml.safe_load(Path(cli_args.config).read_text())


if __name__ == "__main__":
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    cli_args = verify_cli_args(cli_args_parser().parse_args())
    if el := get_callable(cli_args.remaining[0]):
        el(*cli_args.remaining[:-1])

    if ph.global_vars.sim is None:
        raise ValueError(bad_input_error)

    sim = ph.global_vars.sim
    if cli_args.config:
        from copy import deepcopy
        config = read_yaml(cli_args)
        overrides = deepcopy(sim.__dict__)
        overrides.update(config)
        sim.__dict__.update(overrides)

    print_sim_summary(sim)
    live = Simulator(sim, cli_args.monitor)
    live.initialize()
    live.post_advance = PostOp(cli_args, live)
    live.run().reset()
