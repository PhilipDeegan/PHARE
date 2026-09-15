#!/usr/bin/env python3
"""
Shades every level-ghost cell of a level by how many level-ghost particles
land in it:
  white = 0 particles  (a hole -- nobody filled this cell)
  black = max(count) across all level-ghost cells of this level
  grey  = in between   (duplicates pile counts above the expected 1x fill)

Meant to spot duplicate/dropped level-ghost particles at same-level patch/tile
boundaries.

Usage: python3 check_levelghost_duplicates.py <diag_dir> [population] [level]
"""
import sys
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

from pyphare.core.box import grow
from pyphare.pharesee.run import Run
from pyphare.pharesee.geometry import level_ghost_boxes

# how many extra cells of context to show around each tight ghost-box piece
# in the zoomed panels, so the shaded cells aren't squashed flush against the
# panel border with nothing around them to judge scale/position against.
ZOOM_MARGIN = 3

outputpath = Path("phare_outputs/check_levelghost_duplicates")


def per_cell_counts(icells, box):
    counts = np.zeros(box.shape, dtype=int)
    if icells.shape[0] == 0:
        return counts
    for idx in icells - box.lower:
        counts[tuple(idx)] += 1
    return counts


def domain_edges(box, domain):
    """Coordinates where `box` (a ghost-only piece) touches `domain`'s
    boundary, so the plot can draw the domain edge directly instead of
    requiring the viewer to compare axis numbers against the domain box."""
    edges = []
    if box.upper[0] < domain.lower[0]:
        edges.append(("v", domain.lower[0]))
    elif box.lower[0] > domain.upper[0]:
        edges.append(("v", domain.upper[0] + 1))
    if box.upper[1] < domain.lower[1]:
        edges.append(("h", domain.lower[1]))
    elif box.lower[1] > domain.upper[1]:
        edges.append(("h", domain.upper[1] + 1))
    return edges


def side_label(box, domain):
    """Plain-English description of which side/corner of the patch `box`
    is (as opposed to its raw Box(...) repr)."""
    parts = []
    for kind, coord in domain_edges(box, domain):
        if kind == "v":
            parts.append("left" if coord == domain.lower[0] else "right")
        else:
            parts.append("bottom" if coord == domain.lower[1] else "top")
    if not parts:
        return "(not outside domain?)"
    return "-".join(parts) + (" corner" if len(parts) > 1 else " edge")


def main():
    diag_dir = sys.argv[1]
    pop = sys.argv[2] if len(sys.argv) > 2 else "protons"
    ilvl = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    run = Run(diag_dir)
    hier = run.GetParticles(0.0, pop, type="levelGhost")

    gaboxes_per_level = level_ghost_boxes(hier, "levelGhost")
    if ilvl not in gaboxes_per_level:
        raise RuntimeError(f"no level ghost boxes for level {ilvl} in {diag_dir}")
    gaboxes_list = gaboxes_per_level[ilvl][f"{pop}_levelGhost"]

    pdata_key = f"{pop}_levelGhost"
    patch_id_of = {
        id(patch.patch_datas[pdata_key]): patch.id
        for patch in hier.level(ilvl).patches
        if pdata_key in patch.patch_datas
    }

    per_box_counts = []
    for gabox in gaboxes_list:
        pdata = gabox["pdata"]
        for box in gabox["boxes"]:
            icells = pdata.dataset.select(box).iCells
            per_box_counts.append((box, per_cell_counts(icells, box), pdata.box, pdata))

    vmax = max((c.max() for _, c, _, _ in per_box_counts), default=0)
    vmax = max(vmax, 1)  # keep a non-degenerate scale if everything is empty

    fig, ax = plt.subplots()
    im = None
    for box, counts, _, _ in per_box_counts:
        x = np.arange(box.lower[0], box.upper[0] + 2)
        y = np.arange(box.lower[1], box.upper[1] + 2)
        im = ax.pcolormesh(x, y, counts.T, cmap="gray_r", vmin=0, vmax=vmax)

    if im is not None:
        plt.colorbar(im, ax=ax, label="level-ghost particles per cell")
    ax.set_aspect("equal")
    ax.set_title(
        f"L{ilvl} level-ghost particle counts ({pop}), max={vmax}\n"
        "(ghost layer is only a few cells wide -- see _zoomed.png for detail)"
    )

    outputpath.mkdir(parents=True, exist_ok=True)
    outfile = outputpath / f"levelghost_counts_{pop}_L{ilvl}_overview.png"
    fig.savefig(outfile, dpi=200)
    plt.close(fig)
    print(f"saved {outfile} (max count = {vmax})")

    # The ghost layer is only `ghosts` cells wide (often 1), which is a
    # sub-pixel hairline against a whole-domain plot at equal aspect --
    # render each strip on its own axis with auto (non-equal) aspect so it
    # fills the panel and individual cells are actually visible.
    ncols = 4
    nrows = -(-len(per_box_counts) // ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False
    )
    im = None
    for i, (box, counts, domain, pdata) in enumerate(per_box_counts):
        zax = axes[i // ncols][i % ncols]
        # pad the tight ghost-box piece with a few cells of context, clamped
        # to what this patch actually stores (ghost_box) -- the panel would
        # otherwise be flush against the shaded cells with no surrounding
        # data to judge them against.
        padded = grow(box, [ZOOM_MARGIN] * box.ndim) * pdata.ghost_box
        padded_counts = per_cell_counts(pdata.dataset.select(padded).iCells, padded)
        x = np.arange(padded.lower[0], padded.upper[0] + 2)
        y = np.arange(padded.lower[1], padded.upper[1] + 2)
        im = zax.pcolormesh(x, y, padded_counts.T, cmap="gray_r", vmin=0, vmax=vmax)
        for kind, coord in domain_edges(box, domain):
            (zax.axvline if kind == "v" else zax.axhline)(
                coord, color="tab:blue", linewidth=1.5, linestyle="--"
            )

        # small (corner) panels get the exact count written on each cell --
        # color alone can't be read precisely at this size, and this is
        # exactly where the interesting cell-by-cell detail lives.
        if padded.shape[0] <= 10 and padded.shape[1] <= 10:
            for ix in range(padded.shape[0]):
                for iy in range(padded.shape[1]):
                    count = padded_counts[ix, iy]
                    zax.text(
                        padded.lower[0] + ix + 0.5,
                        padded.lower[1] + iy + 0.5,
                        str(count),
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="white" if count > vmax / 2 else "black",
                    )

        patch_id = patch_id_of.get(id(pdata), "?")
        zax.set_title(f"patch {patch_id}: {side_label(box, domain)}", fontsize=8)
    for j in range(len(per_box_counts), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    if im is not None:
        fig.colorbar(im, ax=axes, label="level-ghost particles per cell", shrink=0.6)
    fig.suptitle(
        f"L{ilvl} level-ghost particles per cell, one panel per ghost-region piece ({pop})\n"
        "each patch's ghost border splits into several disconnected pieces (a strip "
        "per side, cells per corner);\n"
        "dashed blue line = where this patch's real domain starts -- numbers/shading "
        "on the near side are the actual ghost-cell counts"
    )

    outfile_zoom = outputpath / f"levelghost_counts_{pop}_L{ilvl}_zoomed.png"
    fig.savefig(outfile_zoom, dpi=150)
    plt.close(fig)
    print(f"saved {outfile_zoom}")


if __name__ == "__main__":
    main()
