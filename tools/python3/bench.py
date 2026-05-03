#
# parsing PHARE scope funtion timers
#

import sys
import ast
import tarfile
import argparse
import numpy as np
from io import BytesIO
from pathlib import Path
from typing import Callable
from dataclasses import dataclass, field

from pyphare.pharesee.run import Run
import phlop.timing.scope_timer as phst


def _as_lines(input):
    if type(input) is bytes:
        return input.decode("utf8").split("\n")
    with open(input, "r") as file:
        return file.readlines()


@dataclass
class LevelStats:
    ilvl: int = 0
    patches: int = 0
    cells: int = 0
    patch_ghost_cells: int = 0
    level_ghost_cells: int = 0
    particles: int = 0
    cell_ratio: int = 0


@dataclass
class SummaryFile:
    kwargs: dict
    levels: list
    particles: int

    def __iter__(self):
        return iter((self.kwargs, self.levels, self.particles))


def _parse_level_info(level):
    trim = [")", ","]
    while any(level.endswith(t) for t in trim):
        level = level[:-1]
    bits = [s for s in level.strip().split(" ") if s]
    ilvl = bits[0].split(":")[1]
    for bit in bits:
        needle = "parts("
        if bit.startswith(needle):
            particles = int(bit[len(needle)])
    bits = bits[1][1:-1].split(",")
    patches, cells, p_ghosts_cells, l_ghosts_cells = [int(s) for s in bits]
    ratio = (p_ghosts_cells + l_ghosts_cells) / cells
    return LevelStats(
        ilvl, patches, cells, p_ghosts_cells, l_ghosts_cells, particles, ratio
    )


def parse_summary_file(input):
    lines = _as_lines(input)
    kwargs = lines[0]
    levels = lines[1:-1]
    total = lines[-1].split(":")[1]
    return SummaryFile(
        ast.literal_eval(kwargs),
        [_parse_level_info(lvl.strip()) for lvl in levels],
        total,
    )


def percent_level_refined(summary_file):
    kwargs, levels, _ = summary_file
    cells = kwargs["cells"]
    ndim = len(cells)
    assert ndim == 3
    dic = {i: 0 for i in range(0, len(levels) - 1)}
    for i in range(0, len(levels) - 1):
        dic[i] = int((levels[i + 1].cells / (2**ndim)) / levels[i].cells * 100)
    return dic


def average_L0_advance_time_from(timer_files):
    advance_times = [tf.advances[0].run_time for tf in timer_files]
    return np.average(advance_times)


def max_memory_required_from(stat_files):
    ret = 0
    for file in stat_files:
        for snapshot in file.snapshots:
            ret = max(ret, snapshot.mem_used_mb)
    return ret


def _print_stats(summary_file, timer_files, stat_files):
    pc_refined = percent_level_refined(summary_file)
    tags = {k: v for k, v in summary_file.kwargs.items() if k.startswith("tag")}
    time = average_L0_advance_time_from(timer_files)
    max_mem = max_memory_required_from(stat_files)
    print(
        f"% refined:{pc_refined}",
        f"N particles:{summary_file.particles}",
        f"cell_ratio:{int(summary_file.levels[1].cell_ratio*100)}",
        f"time:{time / 1e9:.3f}s",
        f"MEM:{max_mem}mb",
        tags,
    )


def extract_nested_tarfile(data):
    with tarfile.open(fileobj=data, mode="r:*") as nestedTar:

        def extract(member):
            if file := nestedTar.extractfile(member):
                return file.read()

        return {
            k: v
            for k, v in {member.name: extract(member) for member in nestedTar}.items()
            if v
        }


def handle_timer_file(data):
    import tools.python3.phloping as phloping

    return phloping.make_scope_timer_file_from(_as_lines(data))


def handle_stat_file(data):
    import tools.python3.stats as stats

    return stats.lines_parser(data)


def handle_pertumation(objects):
    stat_files = [
        handle_stat_file(v) for k, v in objects.items() if "phare/stats/" in k
    ]
    timer_files = [
        handle_timer_file(v) for k, v in objects.items() if "phare/timings/" in k
    ]
    for k, v in objects.items():
        if k.endswith("phare/summary.txt"):
            _print_stats(parse_summary_file(v), timer_files, stat_files)


def extract_nested_directory(path):
    def extract(file):
        with open(file, "rb") as f:
            return f.read()

    return {
        str(d / file): extract(str(d / file))
        for d, _, files in path.walk()
        for file in files
    }


def handle_nested_directory(path):
    handle_pertumation(extract_nested_directory(path))


def handle_nested_tarfile(nested):
    if type(nested) is not BytesIO:
        with open(str(nested), "rb") as f:
            nested = BytesIO(f.read())
    handle_pertumation(extract_nested_tarfile(nested))


def handle_tarfile(path):
    with tarfile.open(str(path), "r:gz") as tar:
        for member in tar:
            if member.name.endswith(".tar.gz"):
                handle_nested_tarfile(BytesIO(tar.extractfile(member).read()))


def this_or_that(pred, a, b, *args, **kwargs):
    return a(*args, **kwargs) if pred else b(*args, **kwargs)


def is_tar_file(any):
    return str(any).endswith(".tar.gz")


def handle_nested_item(path):
    this_or_that(
        is_tar_file(path), handle_nested_tarfile, handle_nested_directory, path
    )


def handle_directory(path):
    for d in Path(path).iterdir():
        handle_nested_item(d)


def handle_input(path):
    this_or_that(is_tar_file(path), handle_tarfile, handle_directory, path)


def print_summary(filepath=None, root_id=None):
    if filepath:
        return handle_input(filepath)

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default=None, help="input")
    args = parser.parse_args()
    if args.input:
        return handle_input(args.input)

    for el in Path(".phare_bench").iterdir():
        ...  # take last
    handle_input(el)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("usage: $function_name -h")
        print(
            "available functions:\n\t"
            + "\n\t".join([k for k, v in globals().items() if k.startswith("print")]),
        )
    elif len(sys.argv) > 1:
        fn = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        if fn not in globals():
            raise ValueError("requested function does not exist")
        try:
            globals()[fn]()
        except BrokenPipeError:
            ...
