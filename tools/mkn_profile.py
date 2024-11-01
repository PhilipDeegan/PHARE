import os
import sys
import subprocess
import numpy as np

from pathlib import Path
import shutil

PATCHES = [1, 2, 3, 4, 5, 10]
CELLS = [5, 6, 7, 8, 9, 10, 12, 15, 20]
PPC = [10, 20, 50, 100]

PATCHES = [10]
CELLS = [5]
PPC = [5, 10]


permutables = [
    ("patches", PATCHES),
    ("cells", CELLS),
    ("ppc", PPC),
]

PHARE_ASYNC_TIMES = ".phare/async/multi_updater"

fn_strings = {
    0: "sync ghosts",
    1: "Boris::move_domain",
    2: "Boris::move_ghosts",
    3: "sync",
    4: "Group_barrier",
    5: "Domain_insert",
    6: "Group_barrier",
    7: "Deposit",
}


def run_permutation(patches, cells, ppc):
    times_dir = f"{PHARE_ASYNC_TIMES}/{patches}/{cells}/{ppc}"
    env = os.environ.copy()
    env_update = {
        "PHARE_ASYNC_TIMES": times_dir,
        "PHARE_PATCHES": f"{patches}",
        "PHARE_CELLS": f"{cells}",
        "PHARE_PPC": f"{ppc}",
    }
    env.update(env_update)
    print({k: v for k, v in env.items() if k.startswith("PHARE_")})

    subprocess.run(
        ["mkn", "run", "-p", "test_core"] + sys.argv[1:], check=True, env=env
    )

    shutil.copy(".phare_times.0.txt", f"{times_dir}/scope_times.txt")
    return times_dir


def parse_times(times_dir, patches, cells, ppc):
    times_per_type = {}

    for f in Path(times_dir).glob("*.txt"):
        type_id = f.stem
        if type_id == "scope_times":
            continue

        times = []
        with open(f) as file:
            for line in file.readlines():
                line = line.strip()
                times.append(line.split(" "))
        times_per_type[type_id] = times

    return {times_dir: (times_per_type, patches, cells, ppc)}


def plot_fn(times, fn):
    import matplotlib.pyplot as plt

    x_axis = [0.1, 0.2, 0.3]

    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    fn_0_times = []

    for times_dir, times_per_type_tuple in times.items():
        times_per_type, patches, cells, ppc = times_per_type_tuple
        for type_id, bits_list in times_per_type.items():
            type_str = type_id.split(",")[1]
            type_str = f"{type_str}/{patches}/{cells}/{ppc}"
            fn_0_times = [[], [], []]

            pid = 0
            for bits in bits_list:
                if bits[1] == f"{fn}":
                    fn_0_times[pid % 3].append(int(bits[2]))
                    pid += 1

            fn_0_times = [np.asarray(ts).mean() for ts in fn_0_times]
            ax.plot(x_axis, fn_0_times, ":", label=f"{type_str}")

    box = ax.get_position()
    ax.set_position(
        [box.x0, box.y0 + box.height * 0.2, box.width * 0.8, box.height * 0.8]
    )

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fancybox=True,
        shadow=True,
    )
    ax.set_title(fn_strings[fn])

    plt.ylabel("nanoseconds")

    ax.set_xticks(x_axis)
    ax.set_xticklabels(
        ["domain", "patchghost", "levelghost"], rotation="vertical", fontsize=18
    )
    fig.savefig(f"profile_plot.{fn}.png")


def plot(times):
    for i in range(len(fn_strings)):
        plot_fn(times, i)


def main():
    import itertools

    times = {}
    permutations = itertools.product(*[e[1] for e in permutables])
    for permutation in permutations:
        times.update(parse_times(run_permutation(*permutation), *permutation))

    plot(times)


if __name__ == "__main__":
    main()
