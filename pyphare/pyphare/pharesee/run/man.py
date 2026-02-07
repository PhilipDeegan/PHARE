#
#
#

from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy import hierarchy_compute as hc

from pyphare.pharesee.hierarchy import patch
from pyphare.pharesee.hierarchy import patchdata

from pyphare.pharesee.hierarchy import ScalarField, VectorField, TensorField

from pyphare.pharesee.hierarchy.tensorfield import AnyTensorField

from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field
from .utils import (
    _compute_to_primal,
    _compute_pop_pressure,
    _compute_pressure,
    _compute_current,
    _compute_divB,
    _get_rank,
    make_interpolator,
)


class LazyFieldData:
    def __init__(self, data, default_mutators, select_mutators):
        self.box = data.box
        self.data = data
        self.patch_datas = {data.name: data}
        self.default_mutators = default_mutators or []
        self.select_mutators = select_mutators or []
        self.loaded = False

    def __getitem__(self, get):
        if not self.loaded:
            get_all = get == slice(None) or get is None
            if get_all:
                for mut in self.default_mutators:
                    self.patch_datas = {ret["name"]: ret["data"] for ret in mut(self)}
            # todo slice
            for mut in self.select_mutators:
                self.patch_datas = {ret["name"]: ret["data"] for ret in mut(self)}
            self.loaded = True
        return self.patch_datas[self.data.name].dataset


def to_lazy_dataset(patch, hinfo, default_mutators, select_mutators):
    return tuple(
        {
            "name": name,
            "data": pd.copy_as(LazyFieldData(pd, default_mutators, select_mutators)),
        }
        for name, pd in patch.patch_datas.items()
    )


class LazyPatch(patch.Patch):
    def __init__(self, patch, mutators):
        super().__init__(patch.patch_datas, patch.patch_id, patch.layout, patch.attrs)
        self.pds = {}
        self.mutators = mutators or []

    def __getitem__(self, key):
        if key not in pds:
            for mut in self.select_mutators:
                self.pds = {ret["name"]: ret["data"] for ret in mut(self)}
        return self.patch_datas[key]


def lazy_patch_hierarchy(hier):
    for time, levels in hier.time_hier.items():
        for ilvl, lvl in levels.items():
            lvl.patches = [LazyPatch(patch) for patch in lvl.patches]
    return hier


class MutatorScalarField(AnyTensorField):
    def __init__(self, hier, default_mutators=None, select_mutators=None):
        super().__init__(
            compute_hier_from(
                to_lazy_dataset,
                hier,
                default_mutators=default_mutators,
                select_mutators=select_mutators,
            )
        )


class MutatorVectorField(AnyTensorField):
    def __init__(self, hier, default_mutators=None, select_mutators=None):
        super().__init__(
            compute_hier_from(
                to_lazy_dataset,
                hier,
                default_mutators=default_mutators,
                select_mutators=select_mutators,
            )
        )


class MutatorTensorField(AnyTensorField):
    def __init__(self, hier, default_mutators=None, select_mutators=None):
        super().__init__(
            compute_hier_from(
                to_lazy_dataset,
                hier,
                default_mutators=default_mutators,
                select_mutators=select_mutators,
            )
        )


class RunFunc:
    def __init__(self, λ, *args, **kwargs):
        self.λ = λ
        self.args = args
        self.kwargs = {**kwargs}

    def __call__(self, *args, **kwargs):
        self.λ(*args, *self.args, **kwargs, **self.kwargs)


class RunMan:
    def __init__(self, run, *args, **kwargs):
        self.run = run
        self.patch_mutators = []
        self.level_mutators = []
        for key, val in kwargs.items():
            setattr(self, key, val)

    def GetVi(self, time, merged=False, interp="nearest"):
        return MutatorVectorField(
            self.run._get_hier_for(time, "ions_bulkVelocity"),
            [hc.drop_ghosts],
            [],
        )

    def GetN(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        return ScalarField(
            MutatorScalarField(
                self.run._get_hier_for(time, f"ions_pop_{pop_name}_density", **kwargs),
                [hc.drop_ghosts],
            )
        )

    def GetFlux(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        return VectorField(
            MutatorVectorField(
                self.run._get_hier_for(time, f"ions_pop_{pop_name}_flux", **kwargs),
                [hc.drop_ghosts],
                [],
            )
        )

    def GetPressure(self, time, pop_name, merged=False, interp="nearest", **kwargs):
        M = MutatorTensorField(
            self.run._get_hierarchy(
                time, f"ions_pop_{pop_name}_momentum_tensor.h5", **kwargs
            ),
            [hc.drop_ghosts],
        )
        V = self.GetFlux(time, pop_name, **kwargs)
        N = self.GetN(time, pop_name, **kwargs)
        return TensorField(
            compute_hier_from(
                _compute_pop_pressure,
                (M, V, N),
                popname=pop_name,
                mass=self.run.GetMass(pop_name, **kwargs),
            )
        )
