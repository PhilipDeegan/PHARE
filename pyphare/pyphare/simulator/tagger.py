"""Tagging outside a simulation: where a run's levels will be, from its initial fields.

    tag_boxes = Tagger(sim).tag()        # {level: [Box]}, level 0 the domain
    weights   = Tagger(sim).weights()    # per level 0 cell, the work expected on it and above
    Tagger(sim).hierarchy().plot(qty="tags", plot_patches=True)  # as a run's hierarchy

No Simulator, no MPI, no particles. Per level, from the coarsest:

  1. the tagged quantities' init functions (the model's bx, by, bz) evaluated over each patch's
     field allocation (ghosts included) at their Yee centring, periodic directions wrapped
  2. the run's own tagging strategy on those fields (cpp_etc.tag_patch_<dim>_<interp>: the
     method, quantities and thresholds of sim's tagging config), per patch
  3. the tagged cells clustered into the next level's boxes as SAMRAI's TileClustering: tagged
     cells grown by tag_buffer, the tag boxes (tag_box_size cells from cell 0, SAMRAI's "tiles":
     the patches to be) holding them, their properly nested part, merged two by two while their
     union is a box, refined

Level 1 is what a run makes at t=0 (its level 0 B is the functions'), given the same patches:
patches are chopped as SAMRAI does on one rank (largest_patch_size, cells shared evenly), and a
criterion whose stencil reaches past the field ghosts leaves the patch borders untagged. Finer
levels sample the functions at their resolution, where a run refines B from the coarser level:
an estimate there, close but not identical near the threshold.

Only quantities known before a run can be tagged here: B (hybrid). Particle moments, E and J
are not (pharein rejects them for tagging anyway).
"""

from dataclasses import dataclass, field

import numpy as np

from pyphare.core.box import Box
from pyphare.core.gridlayout import HybridGridLayoutFor

_TAGGING_PATH = "simulation/AMR/refinement/tagging/"


class _DictCollector:
    """a DictPopulator collecting into a nested dict what a run writes under _TAGGING_PATH"""

    def __init__(self):
        self.dict = {}

    def _set(self, path, value):
        if not path.startswith(_TAGGING_PATH):
            return
        keys = path[len(_TAGGING_PATH) :].strip("/").split("/")
        d = self.dict
        for k in keys[:-1]:
            d = d.setdefault(k, {})
        d[keys[-1]] = value

    def add_int(self, path, val):
        self._set(path, int(val))

    def add_double(self, path, val):
        self._set(path, float(val))

    def add_string(self, path, val):
        self._set(path, str(val))

    def add_bool(self, path, val):
        self._set(path, bool(val))


def tagging_dict(sim):
    """the dict a run reads at simulation/AMR/refinement/tagging"""
    if sim.refinement != "tagging":
        raise ValueError("Tagger: the simulation does not refine by tagging")
    dp = _DictCollector()
    dp.add_double(_TAGGING_PATH + "threshold", sim.tagging_threshold)
    dp.add_string(_TAGGING_PATH + "hybrid_method", sim.tagging.method)
    sim.tagging.populate_dict(dp)
    return dp.dict


def chop_evenly(box, size, granularity=1):
    """SAMRAI's patches of a box on one rank: as many as size requires, the cells shared evenly
    in units of granularity cells, the leftover units one each to the first patches"""
    cuts = []
    for d in range(len(box.lower)):
        width = int(box.upper[d] - box.lower[d] + 1)
        n = -(-width // int(size[d]))
        units = -(-width // granularity)
        cuts.append(
            [int(box.lower[d]) + granularity * (-(-k * units // n)) for k in range(n)]
            + [int(box.upper[d]) + 1]
        )
    pieces = [([], [])]
    for c in cuts:
        pieces = [
            (lo + [c[k]], up + [c[k + 1] - 1])
            for lo, up in pieces
            for k in range(len(c) - 1)
        ]
    return [Box(np.array(lo), np.array(up)) for lo, up in pieces]


def sweep_boxes(mask):
    """the True cells of a bool array as disjoint boxes [(lower, upper)]: runs along the first
    direction, then equal runs merged along the next directions in turn"""
    ndim = mask.ndim
    boxes = []  # (lower, upper) arrays
    idx = np.argwhere(mask)
    if not len(idx):
        return boxes
    # runs along direction 0, per line of the other directions
    order = np.lexsort([idx[:, 0]] + [idx[:, d] for d in range(1, ndim)])
    idx = idx[order]
    runs = []
    start = idx[0].copy()
    prev = idx[0].copy()
    for c in idx[1:]:
        if np.array_equal(c[1:], prev[1:]) and c[0] == prev[0] + 1:
            prev = c
            continue
        runs.append((start.copy(), prev.copy()))
        start, prev = c.copy(), c.copy()
    runs.append((start, prev))
    # merge along directions 1.., boxes equal in every other direction and adjacent in d
    for d in range(1, ndim):
        others = [i for i in range(ndim) if i != d]
        runs.sort(
            key=lambda b: (
                tuple(b[0][i] for i in others)
                + tuple(b[1][i] for i in others)
                + (b[0][d],)
            )
        )
        merged = []
        for lo, up in runs:
            if merged:
                mlo, mup = merged[-1]
                same = all(mlo[i] == lo[i] and mup[i] == up[i] for i in others)
                if same and mup[d] + 1 == lo[d]:
                    mup[d] = up[d]
                    continue
            merged.append((lo.copy(), up.copy()))
        runs = merged
    boxes = runs
    return boxes


def merge_boxes(boxes):
    """disjoint boxes [(lower, upper)] merged while two are adjacent along a direction with the
    same extent in every other: passes over the directions until none merge"""
    boxes = [(np.asarray(lo).copy(), np.asarray(up).copy()) for lo, up in boxes]
    if not boxes:
        return boxes
    ndim = len(boxes[0][0])
    while True:
        before = len(boxes)
        for d in range(ndim):
            others = [i for i in range(ndim) if i != d]
            boxes.sort(
                key=lambda b: (
                    tuple(b[0][i] for i in others)
                    + tuple(b[1][i] for i in others)
                    + (b[0][d],)
                )
            )
            merged = []
            for lo, up in boxes:
                if merged:
                    mlo, mup = merged[-1]
                    if (
                        all(mlo[i] == lo[i] and mup[i] == up[i] for i in others)
                        and mup[d] + 1 == lo[d]
                    ):
                        mup[d] = up[d]
                        continue
                merged.append((lo, up))
            boxes = merged
        if len(boxes) == before:
            return boxes


@dataclass
class TaggedLevel:
    number: int
    boxes: list  # [Box], this level's index space
    tagged: int = 0  # cells tagged on it (they made the next level)

    @property
    def cells(self):
        return int(sum(b.nCells() for b in self.boxes))


@dataclass
class Tagger:
    sim: object
    tag_box_size: int = 8  # SAMRAI TileClustering's default tile size
    levels: list = field(default_factory=list)
    tags: dict = field(default_factory=dict)  # {level: [(patch Box, tagged cells)]}
    keep_tags: bool = True  # the patches' tags, for hierarchy() (off for large domains)

    def __post_init__(self):
        sim = self.sim
        self.ndim = sim.ndim
        nd = (self.ndim,)
        self.cells = np.broadcast_to(np.asarray(sim.cells, dtype=int), nd).copy()
        self.dl = np.broadcast_to(np.asarray(sim.dl, dtype=float), nd).copy()
        self.origin = np.broadcast_to(
            np.asarray(sim.origin or 0.0, dtype=float), nd
        ).copy()
        boundary = getattr(sim, "boundary_types", None) or "periodic"
        if isinstance(boundary, str):
            boundary = [boundary] * self.ndim
        self.periodic = np.asarray([b == "periodic" for b in boundary])
        self.ratio = int(getattr(sim, "refinement_ratio", 2))
        self.interp = sim.interp_order
        self.dict = tagging_dict(sim)
        model = sim.model.model_dict
        self.functions = {"Bx": model["bx"], "By": model["by"], "Bz": model["bz"]}
        size = getattr(sim, "largest_patch_size", None)
        self.patch_size = np.broadcast_to(
            np.asarray(size if size is not None else self.cells * 2**30, dtype=int), nd
        ).copy()
        self.tag_buffer = int(getattr(sim, "tag_buffer", 1))
        self.nesting_buffer = int(np.max(getattr(sim, "nesting_buffer", 1)))

    def _cpp(self):
        from pyphare import cpp

        fn = getattr(cpp.cpp_etc_lib(), f"tag_patch_{self.ndim}_{self.interp}", None)
        if fn is None:
            raise RuntimeError("Tagger: cpp_etc has no tag_patch (an older build?)")
        return fn

    def domain(self, ln):
        return Box(np.zeros(self.ndim, dtype=int), self.cells * self.ratio**ln - 1)

    def _values(self, name, ln, patch):
        """the component over its allocation on patch (the hybrid layout's coordinates, ghosts
        included), periodic directions wrapped into the domain"""
        dx = self.dl / self.ratio**ln
        origin = self.origin + patch.lower * dx
        layout = HybridGridLayoutFor(patch, origin.tolist(), dx.tolist(), self.interp)
        length = self.dl * self.cells
        xs = []
        for d, direction in enumerate("xyz"[: self.ndim]):
            x = np.asarray(layout.yeeCoordsFor(name, direction, withGhosts=True))
            if self.periodic[d]:
                x = self.origin[d] + np.mod(x - self.origin[d], length[d])
            xs.append(x)
        mesh = np.meshgrid(*xs, indexing="ij")
        v = np.asarray(self.functions[name](*[m.ravel() for m in mesh]), dtype=float)
        return np.array(
            np.broadcast_to(v, (mesh[0].size,)).reshape(mesh[0].shape), order="C"
        )

    def _tags(self, ln, patch):
        dx = self.dl / self.ratio**ln
        b = [self._values(n, ln, patch) for n in ("Bx", "By", "Bz")]
        tags = self._cpp()(
            self.dict, dx.tolist(), patch.lower.tolist(), patch.upper.tolist(), *b
        )
        return np.asarray(tags).reshape(patch.shape, order="F").astype(bool)

    def _patches(self, ln, box):
        size, gran = self.patch_size, 1
        if ln > 0:  # refined: patches of whole coarser cells
            gran = self.ratio
            size = np.maximum(gran, self.patch_size // gran * gran)
        return chop_evenly(box, size, gran)

    def _grid_shape(self, ln):
        return -(-np.asarray(self.domain(ln).shape) // self.tag_box_size)

    def _mark(self, grown, ln, patch, t):
        """the tag boxes (level ln's grid) the patch's tagged cells grown by tag_buffer touch"""
        shape = np.asarray(self.domain(ln).shape)
        size, tb = self.tag_box_size, self.tag_buffer
        cells = np.argwhere(t) + patch.lower
        if not len(cells):
            return
        sides = []  # per direction: the tag boxes of the lowest and highest grown cell
        for d in range(self.ndim):
            lo, hi = cells[:, d] - tb, cells[:, d] + tb
            if self.periodic[d]:
                lo, hi = np.mod(lo, shape[d]), np.mod(hi, shape[d])
            else:
                lo, hi = np.clip(lo, 0, shape[d] - 1), np.clip(hi, 0, shape[d] - 1)
            sides.append(
                (lo // size, hi // size)
            )  # tag_buffer < size: at most these two
        for pick in np.ndindex(*([2] * self.ndim)):
            grown[tuple(sides[d][pick[d]] for d in range(self.ndim))] = True

    def _coverage(self, ln, level_boxes):
        """per tag box of level ln, its cells on the level's boxes (boxes do not overlap), and
        its cells in the domain"""
        shape = np.asarray(self.domain(ln).shape)
        grid = self._grid_shape(ln)
        size = self.tag_box_size

        def per_axis(lo, up, d):  # cells of [lo, up] in each tag box along d
            k = np.arange(grid[d])
            return np.maximum(
                0,
                np.minimum(up, np.minimum(k * size + size - 1, shape[d] - 1))
                - np.maximum(lo, k * size)
                + 1,
            )

        def outer(counts):
            out = counts[0]
            for c in counts[1:]:
                out = np.multiply.outer(out, c)
            return out

        covered = np.zeros(grid, dtype=np.int64)
        for b in level_boxes:
            covered += outer(
                [per_axis(b.lower[d], b.upper[d], d) for d in range(self.ndim)]
            )
        cells = outer([per_axis(0, shape[d] - 1, d) for d in range(self.ndim)])
        return covered, cells

    def _nested_cells(self, ln, level_boxes, lower, upper):
        """the cells lower..upper whose every cell within nesting_buffer is on level ln's boxes
        (periodic images wrapped, past a non periodic boundary counting as covered): a mask"""
        shape = np.asarray(self.domain(ln).shape)
        nb = self.nesting_buffer
        glo, gup = lower - nb, upper + nb
        covered = np.zeros(gup - glo + 1, dtype=bool)
        for d in range(self.ndim):  # outside a non periodic domain: covered
            if not self.periodic[d]:
                idx = np.arange(glo[d], gup[d] + 1)
                out = (idx < 0) | (idx >= shape[d])
                view = [None] * self.ndim
                view[d] = slice(None)
                covered |= out[tuple(view)]
        shifts = [np.asarray(s) * shape for s in np.ndindex(*([3] * self.ndim))]
        for b in level_boxes:
            for s in shifts:
                s = np.where(self.periodic, s - shape, 0)
                lo = np.maximum(glo, b.lower + s)
                up = np.minimum(gup, b.upper + s)
                if np.all(lo <= up):
                    covered[
                        tuple(slice(l - g, u - g + 1) for l, u, g in zip(lo, up, glo))
                    ] = True
        nested = covered.copy()
        for d in range(self.ndim):
            out = nested.copy()
            for k in range(1, nb + 1):
                out[
                    tuple(
                        slice(k, None) if i == d else slice(None)
                        for i in range(self.ndim)
                    )
                ] &= nested[
                    tuple(
                        slice(None, -k) if i == d else slice(None)
                        for i in range(self.ndim)
                    )
                ]
                out[
                    tuple(
                        slice(None, -k) if i == d else slice(None)
                        for i in range(self.ndim)
                    )
                ] &= nested[
                    tuple(
                        slice(k, None) if i == d else slice(None)
                        for i in range(self.ndim)
                    )
                ]
            nested = out
        return nested[
            tuple(slice(nb, -nb) if nb else slice(None) for _ in range(self.ndim))
        ]

    def _cluster(self, ln, level_boxes, grown):
        """level ln + 1's boxes from level ln's grown tags on the tag box grid"""
        shape = np.asarray(self.domain(ln).shape)
        size = self.tag_box_size
        covered, cells = self._coverage(ln, level_boxes)
        full = covered == cells

        # wholly nested: the tag box and every neighbour within one tag box (nesting_buffer <
        #  tag_box_size) fully covered (past a non periodic boundary: covered)
        whole = full.copy()
        for d in range(self.ndim):
            out = whole.copy()
            for sign in (-1, 1):
                shifted = np.roll(whole, -sign, axis=d)
                if not self.periodic[d]:
                    idx = [slice(None)] * self.ndim
                    idx[d] = slice(-1, None) if sign > 0 else slice(0, 1)
                    shifted[tuple(idx)] = True
                out &= shifted
            whole = out

        boxes = [
            Box(lo * size, np.minimum((up + 1) * size, shape) - 1)
            for lo, up in sweep_boxes(grown & whole)
        ]
        for t in np.argwhere(grown & ~whole & (covered > 0)):  # at the level's edge
            lower = t * size
            upper = np.minimum(lower + size, shape) - 1
            nested = self._nested_cells(ln, level_boxes, lower, upper)
            boxes += [Box(lower + lo, lower + up) for lo, up in sweep_boxes(nested)]
        r = self.ratio
        return [
            Box(lo * r, (up + 1) * r - 1)
            for lo, up in merge_boxes((b.lower, b.upper) for b in boxes)
        ]

    def tag(self, depth=None):
        """{level: [Box]}: the levels from the initial fields, up to depth levels above 0
        (default max_nbr_levels - 1)"""
        max_levels = self.sim.max_nbr_levels
        if depth is not None:
            max_levels = min(1 + depth, max_levels)
        self.levels = [TaggedLevel(0, [self.domain(0)])]
        self.tags = {}
        while len(self.levels) < max_levels:
            lvl = self.levels[-1]
            ln = lvl.number
            grown = np.zeros(self._grid_shape(ln), dtype=bool)
            for box in lvl.boxes:
                for patch in self._patches(ln, box):
                    t = self._tags(ln, patch)
                    lvl.tagged += int(t.sum())
                    self._mark(grown, ln, patch, t)
                    if self.keep_tags:
                        self.tags.setdefault(ln, []).append((patch, t))
            if not lvl.tagged:
                break
            boxes = self._cluster(ln, lvl.boxes, grown)
            if not boxes:
                break
            self.levels.append(TaggedLevel(lvl.number + 1, boxes))
        return {lvl.number: lvl.boxes for lvl in self.levels}

    def hierarchy(self):
        """the levels as a pharesee PatchHierarchy (time 0): each level's patches as the tagging
        chopped them, a cell centred "tags" field on each (1 on tagged cells, 0 elsewhere; NaN on
        the finest level, not tagged: transparent in plots), for pharesee's plots"""
        from pyphare.core.gridlayout import GridLayout
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy
        from pyphare.pharesee.hierarchy.patch import Patch
        from pyphare.pharesee.hierarchy.patchdata import FieldData
        from pyphare.pharesee.hierarchy.patchlevel import PatchLevel

        if not self.levels:
            self.tag()
        centering = ["dual"] * self.ndim
        levels = {}
        for lvl in self.levels:
            ln = lvl.number
            dx = self.dl / self.ratio**ln
            tagged = self.tags.get(ln) or [
                (p, np.full(p.shape, np.nan))
                for box in lvl.boxes
                for p in self._patches(ln, box)
            ]
            patches = []
            for i, (box, t) in enumerate(tagged):
                layout = GridLayout(
                    box,
                    (self.origin + box.lower * dx).tolist(),
                    dx.tolist(),
                    self.interp,
                )
                data = FieldData(
                    layout,
                    "tags",
                    t.astype(float),
                    centering=centering,
                    ghosts_nbr=[0] * self.ndim,
                )
                patches.append(Patch({"tags": data}, patch_id=f"p{ln}#{i}"))
            levels[ln] = PatchLevel(ln, patches)
        return PatchHierarchy({0.0: levels}, self.domain(0))

    def weights(self, substeps=True, block=1):
        """per level 0 block (block cells a side; 1: per cell) the work expected on it and above:
        every level's cells over it times the particles per cell, times the level's substeps
        (ratio^l) if substeps"""
        if not self.levels:
            self.tag()
        model = self.sim.model
        ppc = (
            sum(model.model_dict[p]["nbrParticlesPerCell"] for p in model.populations)
            or 1
        )
        blocks = -(-self.cells // block)
        w = np.zeros(blocks)
        for lvl in self.levels:
            r = self.ratio**lvl.number
            scale = ppc * (r if substeps else 1)
            for b in lvl.boxes:
                counts = [
                    np.bincount(np.arange(lo, up + 1) // (r * block), minlength=n)[:n]
                    for lo, up, n in zip(b.lower, b.upper, blocks)
                ]
                per_cell = counts[0]
                for c in counts[1:]:
                    per_cell = np.multiply.outer(per_cell, c)
                w += scale * per_cell
        return w
