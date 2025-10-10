#
#
#

import json
import importlib
from pyphare.cpp import validate

__all__ = ["validate"]

_libs = {}


def simulator_id(sim):
    return f"{sim.ndim}_{sim.interp_order}_{sim.refined_particle_nbr}"


def cpp_lib(sim):
    global _libs
    sim_str = simulator_id(sim)
    if sim_str in _libs:
        return _libs[sim_str]

    mod_str = f"pybindlibs.cpp_{sim_str}"
    _libs[mod_str] = importlib.import_module(mod_str)
    return _libs[mod_str]


def cpp_etc_lib():
    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(sim):
    return getattr(cpp_lib(sim), "Splitter")


def create_splitter(dim, interp, nbrRefinedPart):
    return splitter_type(dim, interp, nbrRefinedPart)()


def split_pyarrays_fn(sim):
    return getattr(cpp_lib(sim), "split_pyarray_particles")


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()
