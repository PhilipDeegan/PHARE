#
#
#


from . import compute as hc
from . import hierarchy_utils as hootils
from .hierarchy import HierarchyAccessor, PatchHierarchy


class AnyTensorField(PatchHierarchy):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def FROM(klass, hier):
        return klass(
            hier.time_hier,
            hier.domain_box,
            hier.refinement_ratio,
            hier.data_files,
            hier.selection_box,
            hier.ephemerals,
        )

    def __getitem__(self, input):
        if type(input) is HierarchyAccessor:
            return super().__getitem__(input)
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)


class TensorField(AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = ["xx", "xy", "xz", "yy", "yz", "zz"]

    @classmethod
    def FROM(cls, hier):
        return super().FROM(cls, hier)

    def __mul__(self, other):
        if type(other) is TensorField:
            raise ValueError(
                "TensorField * TensorField is ambiguous, use pyphare.pharesee.hierarchy.compute.dot or .prod"
            )
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_mul, self, operand=other)
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_add, self, operand=other)
        )

    def __sub__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_sub, self, operand=other)
        )

    def __truediv__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_truediv, self, operand=other)
        )
