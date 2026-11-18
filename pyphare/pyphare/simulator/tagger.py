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

_TAGGING_PATH = "simulation/AMR/refinement/tagging/"

# hybrid Yee layout: per B component, primal or not per direction
_yee_primal = {
    "Bx": (True, False, False),
    "By": (False, True, False),
    "Bz": (False, False, True),
}

# hybrid field ghost width per interpolation order (HybridFieldOptions)
_field_ghost_width = {1: 2, 2: 4, 3: 4}


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


def _union(a, b):
    differ = [
        d
        for d in range(len(a.lower))
        if a.lower[d] != b.lower[d] or a.upper[d] != b.upper[d]
    ]
    if not differ:
        return a
    if len(differ) > 1:
        return None
    d = differ[0]
    if a.upper[d] + 1 != b.lower[d] and b.upper[d] + 1 != a.lower[d]:
        return None
    return Box(np.minimum(a.lower, b.lower), np.maximum(a.upper, b.upper))


def coalesce(boxes):
    """boxes merged two by two while their union is a box, in lower corner order"""
    boxes = list(boxes)
    merged = True
    while merged:
        merged = False
        boxes.sort(key=lambda b: (tuple(b.lower), tuple(b.upper)))
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                u = _union(boxes[i], boxes[j])
                if u is not None:
                    boxes[i] = u
                    del boxes[j]
                    merged = True
                    break
            if merged:
                break
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
        self.ghosts = _field_ghost_width[self.interp]
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
        """the component over its allocation on patch: shape the layout's allocSize"""
        primal = _yee_primal[name][: self.ndim]
        dx = self.dl / self.ratio**ln
        length = self.dl * self.cells
        xs = []
        for d in range(self.ndim):
            idx = np.arange(
                patch.lower[d] - self.ghosts,
                patch.upper[d] + self.ghosts + (2 if primal[d] else 1),
                dtype=float,
            )
            x = (idx + (0.0 if primal[d] else 0.5)) * dx[d]
            if self.periodic[d]:
                x = np.mod(x, length[d])
            xs.append(self.origin[d] + x)
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

    def _cluster(self, ln, level_boxes, tags):
        shape = np.asarray(self.domain(ln).shape)
        size = np.full(self.ndim, self.tag_box_size)

        grown = np.zeros(shape, dtype=bool)
        tb = self.tag_buffer
        for patch, t in tags:
            cells = np.argwhere(t) + patch.lower
            for off in np.ndindex(*([2 * tb + 1] * self.ndim)):
                c = cells + np.asarray(off) - tb
                keep = np.ones(len(c), dtype=bool)
                for d in range(self.ndim):
                    if self.periodic[d]:
                        c[:, d] = np.mod(c[:, d], shape[d])
                    else:
                        keep &= (c[:, d] >= 0) & (c[:, d] < shape[d])
                grown[tuple(c[keep].T)] = True

        covered = np.zeros(shape, dtype=bool)
        for b in level_boxes:
            covered[tuple(slice(lo, up + 1) for lo, up in zip(b.lower, b.upper))] = True
        nested = covered.copy()
        for d in range(self.ndim):
            out = nested.copy()
            for k in range(1, self.nesting_buffer + 1):
                for sign in (-1, 1):
                    shifted = np.roll(nested, -sign * k, axis=d)
                    if not self.periodic[d]:
                        idx = [slice(None)] * self.ndim
                        idx[d] = slice(-k, None) if sign > 0 else slice(0, k)
                        shifted[tuple(idx)] = True
                    out &= shifted
            nested = out

        boxes = []
        for tag_box in sorted({tuple(t) for t in np.argwhere(grown) // size}):
            lower = np.asarray(tag_box) * size
            upper = np.minimum(lower + size - 1, shape - 1)
            sub = nested[tuple(slice(lo, up + 1) for lo, up in zip(lower, upper))]
            if sub.all():
                boxes.append(Box(lower, upper))
            else:
                boxes += [Box(c.copy(), c.copy()) for c in np.argwhere(sub) + lower]
        r = self.ratio
        return [Box(b.lower * r, (b.upper + 1) * r - 1) for b in coalesce(boxes)]

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
            tags = [
                (patch, self._tags(lvl.number, patch))
                for box in lvl.boxes
                for patch in self._patches(lvl.number, box)
            ]
            self.tags[lvl.number] = tags
            lvl.tagged = int(sum(t.sum() for _, t in tags))
            if not lvl.tagged:
                break
            boxes = self._cluster(lvl.number, lvl.boxes, tags)
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

    def weights(self, substeps=True):
        """per level 0 cell the work expected on it and above: every level's cells over it times
        the particles per cell, times the level's substeps (ratio^l) if substeps"""
        if not self.levels:
            self.tag()
        model = self.sim.model
        ppc = (
            sum(model.model_dict[p]["nbrParticlesPerCell"] for p in model.populations)
            or 1
        )
        w = np.zeros(self.cells)
        for lvl in self.levels:
            r = self.ratio**lvl.number
            scale = ppc * (r if substeps else 1)
            for b in lvl.boxes:
                counts = [
                    np.bincount(np.arange(lo, up + 1) // r, minlength=n)[:n]
                    for lo, up, n in zip(b.lower, b.upper, self.cells)
                ]
                per_cell = counts[0]
                for c in counts[1:]:
                    per_cell = np.multiply.outer(per_cell, c)
                w += scale * per_cell
        return w
