#
#
#

import os
import yaml
import inspect
import importlib
import numpy as np
from time import sleep
from pathlib import Path
from dataclasses import dataclass

from pyphare import cpp


USAGE = """!PHARE CLI!


example: python3 pyphare/pyphare --watch 1 tests.functional.harris.harris_2d

Executes Simulation object defined in config() function in file tests/functional/harris/harris_2d.py
            with realtime plots showing live data

"""


@dataclass
class CliHelp:
    input = "Simple YAML overrides for root simulation parameters"
    config = "Simulation configuration function in selected module"
    watch = "Interval to watch realtime plots, default 0 is disabled"
    monitor = "Interval to record system utility information per rank"
    update = "Optional post advance function for live data modifications"
    tui = "Use tui!"


def cli_args_parser():
    import argparse

    _help = CliHelp()
    parser = argparse.ArgumentParser(
        description=USAGE, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-t", "--tui", action="store_true", default=False, help=_help.tui
    )
    parser.add_argument("--input", default="", help=_help.input)
    parser.add_argument("--config", default="config", help=_help.config)
    parser.add_argument("--watch", default=0, help=_help.watch)
    parser.add_argument("--monitor", default=0, help=_help.monitor)
    parser.add_argument("--update", default=None, help=_help.update)
    parser.add_argument("remaining", nargs=argparse.REMAINDER)
    return parser


def verify_cli_args(cli_args):
    try:
        cli_args.watch = int(cli_args.watch)
        cli_args.monitor = int(cli_args.monitor)
    except ValueError:
        raise ValueError("Interval must be an integer")
    return cli_args


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


def get_module(cli_args):
    remaining = cli_args.remaining[0]
    if remaining.endswith(".py"):
        remaining = remaining[:-3]
        remaining = remaining.replace(os.path.sep, ".")
    remaining = f"{remaining}.{cli_args.config}"
    print("remaining", remaining)
    if is_module(remaining):
        return remaining


def get_callable(cli_args):
    path = cli_args.remaining[0]
    if path.endswith(".py"):
        path = path[:-3]
        path = path.replace(os.path.sep, ".")
    if not is_module(path, rethrow=True):
        raise ValueError(bad_input_error)
    for name, el in inspect.getmembers(importlib.import_module(path)):
        if name == cli_args.config:
            return el
    raise ValueError(bad_input_error)


def find_update_in(cli_args, callable):
    for name, el in inspect.getmembers(inspect.getmodule(callable)):
        if name == cli_args.update:
            return el


def realtime_plots(new_time, post_op):
    _locals = dict(idx=0)

    def add_or_get(i):
        if len(post_op.axes) <= i:
            post_op.axes.append(
                post_op.fig.add_subplot(post_op.rows, post_op.cols, i + 1)
            )
        return post_op.axes[i]

    def _try_plot_axis(qty):
        try:
            ax = add_or_get(_locals["idx"])
            hierarchy_from_sim(post_op.live, qty=qty).plot(
                ax=ax, qty=qty, plot_patches=True
            )
            _locals["idx"] += 1
        except FileNotFoundError:
            ...

    def _try_plot_vecfield(qty):
        for c in ["x", "y", "z"]:
            _try_plot_axis(f"{qty}_{c}")

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
    def __init__(self, cli_args, live, update=None):
        import matplotlib.pyplot as plt

        if cli_args.watch:
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
        self.update = update

    def __call__(self, new_time):
        if self.update:
            self.update(self)
        if cli_args.watch and self.tidx % cli_args.watch == 0:
            realtime_plots(new_time, self)
        self.tidx += 1


def humanize_bytes(nbytes):
    suffixes = ["B", "KB", "MB", "GB", "TB", "PB"]
    i = 0
    while nbytes >= 1024 and i < len(suffixes) - 1:
        nbytes /= 1024.0
        i += 1
    f = ("%.2f" % nbytes).rstrip("0").rstrip(".")
    return "%s %s" % (f, suffixes[i])


def memory_required(sim):
    # assume hybrid L0

    particle_size = [56, 64, 76][len(sim.cells) - 1]
    ppc = sum(
        sim.model.model_dict[pop]["nbrParticlesPerCell"]
        for pop in sim.model.populations
    )
    return humanize_bytes(np.prod(sim.cells) * ppc * particle_size)


def print_sim_summary(sim):
    if cpp.mpi_rank() > 0:
        return

    print("Simulation summary:")
    print(f" cells/dl: {sim.cells} / {sim.dl}")
    print(f" final_time: {sim.final_time}")
    print(f" time_step: {sim.time_step}")
    print(f" time_step_nbr: {sim.time_step_nbr}")
    print(f" lower bound memory requirement: {memory_required(sim)}")


import sys

from textual.app import App, ComposeResult
from textual.widgets import Footer, Header, Log
from textual.containers import Container, Horizontal, VerticalScroll
from textual.widgets import Placeholder
from textual.widgets import Input, Static

from textual import work
from textual.worker import Worker
from threading import Thread


class TextualApp(App):
    """An app with a simple log."""

    def __init__(self, max_lines=30):
        super().__init__()
        self.max_lines = max_lines

        # self.simulator = simulator
        self.std_log = Log(id="std_log", max_lines=max_lines * 2, auto_scroll=True)

    def compose(self) -> ComposeResult:
        yield Header()
        yield VerticalScroll(
            Container(
                Placeholder("This is a custom label for p1.", id="p1"),
                Placeholder("Placeholder p2 here!", id="p2"),
                Placeholder(id="p3"),
                Placeholder(id="p4"),
                Placeholder(id="p5"),
                Placeholder(),
                Horizontal(
                    Placeholder(variant="size", id="col1"),
                    Placeholder(variant="text", id="col2"),
                    Placeholder(variant="size", id="col3"),
                    id="c1",
                ),
                id="bot",
            ),
            Container(
                Placeholder(variant="text", id="left"),
                Placeholder(variant="size", id="topright"),
                Placeholder(variant="text", id="botright"),
                id="top",
            ),
            id="content",
        )
        yield self.std_log
        yield Footer()

    def print(self, *a, **kw):
        """redirects to the UI's default std output log"""
        if self.std_log.line_count > self.max_lines:
            self.std_log.refresh_lines(self.max_lines, self.max_lines * 2)
        print(*a, file=self, **kw)

    def flush(self):
        pass

    def write(self, s):
        """for redirecting stdout to a file-like"""
        if self.is_running:
            return self.call_from_anywhere(
                self.std_log.write, s
            )  # std_log is a TextLog
        else:
            return sys.__stdout__.write(s)

    def call_from_anywhere(self, f, *a, **kw):
        try:
            return self.call_from_thread(f, *a, **kw)
        except RuntimeError as e:
            return f(*a, **kw)

    # def watch_data(self) -> None:
    #     sleep(1)
    #     # raise "watch_data"
    #     # self.query_one(Label).update(f"Data: {self.data}")

    # def on_input_submitted(self, message):
    #     sleep(1)
    #     # raise "on_input_submitted"
    #     # self.data = message.value

    # @work(exclusive=True, thread=True)
    # def run(self) -> None:
    #     self.simulator.run()

    # async def on_input_changed(self, message: Input.Changed) -> None:
    #     """Called when the input changes"""
    #     await self.run()


def read_yaml(cli_args):
    return yaml.safe_load(Path(cli_args.input).read_text())


if __name__ == "__main__":
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator, startMPI
    from pyphare.simulator import simulator

    startMPI()

    cli_args = verify_cli_args(cli_args_parser().parse_args())

    update = None
    if el := get_callable(cli_args):
        el()
        if cli_args.update:
            update = find_update_in(cli_args, el)

    if ph.global_vars.sim is None:
        raise ValueError(bad_input_error)

    sim = ph.global_vars.sim
    if cli_args.input:
        from copy import deepcopy

        yaml = read_yaml(cli_args)
        overrides = deepcopy(sim.__dict__)
        overrides.update(yaml)
        sim.__dict__.update(overrides)

    print_sim_summary(sim)
    live = Simulator(sim, print_one_line=True)
    live.initialize()
    live.post_advance = PostOp(cli_args, live, update)

    if cli_args.tui and cpp.mpi_rank() == 0:
        tui = TextualApp()
        simulator.print = tui.print

        t = Thread(target=live.run, daemon=True, args=[cli_args.monitor])
        t.start()
        # sleep(1)
        tui.run()
        # t.stop()
    else:
        # f = open(os.devnull, "w")
        # sys.stdout = f
        # sys.stderr = f
        # sys.stdin = f
        live.run(cli_args.monitor).reset()
