#
#
#

import os
import json
import importlib
from . import validate

__all__ = ["validate"]

_libs = {}


def simulator_id(sim):
    layout = sim.particle_layout
    allocator = sim.allocator

    if getattr(sim, "mhd_timestepper", None):
        Hall = "true" if sim.hall else "false"
        Res = "true" if sim.res else "false"
        Hyper_Res = "true" if sim.hyper_res else "false"
        return (
            f"{sim.ndim}_{sim.interp_order}_{sim.refined_particle_nbr}_{layout}_{allocator}_"
            f"{sim.mhd_timestepper}_{sim.reconstruction}_{sim.limiter}_"
            f"{sim.riemann}_{Hall}_{Res}_{Hyper_Res}"
        )

    return "_".join(
        str(s)
        for s in [
            sim.ndim,
            sim.interp_order,
            sim.refined_particle_nbr,
            layout,
            allocator,
        ]
    )


def cpp_lib(sim):
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim)}"
    if mod_str not in _libs:
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


def split_pyarrays_fn(sim):
    return getattr(cpp_lib(sim), "split_pyarray_particles")


def mpi_initialized():
    return getattr(cpp_etc_lib(), "mpi_initialized")()


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()


def print_rank0(*args, **kwargs):
    def should_print():
        try:
            if mpi_initialized():
                return mpi_rank() == 0
        except ImportError:
            # missing module or mpi not initialized
            ...
        envs = ["OMPI_COMM_WORLD_RANK", "SLURM_PROCID"]
        for env in envs:
            if env in os.environ:
                return int(os.environ[env]) == 0
        return True  # FALL BACK ALWAYS PRINT

    if should_print():
        print(*args, **kwargs)


def layout_mode_enum():
    return getattr(cpp_etc_lib(), "LayoutMode")()


def supported_particle_layouts():
    return getattr(cpp_etc_lib(), "supported_layouts")()


def simulation_to_simopts(sim):
    lib = cpp_etc_lib()
    opts = lib.SimOpts()
    opts.dimension = sim.ndim
    opts.interp_order = sim.interp_order
    opts.nbRefinedPart = sim.refined_particle_nbr
    opts.layout_mode = getattr(lib.LayoutMode, sim.particle_layout)
    opts.alloc_mode = getattr(lib.AllocatorMode, sim.allocator)

    if getattr(sim, "mhd_timestepper", None):
        opts.time_integrator_type = getattr(lib.TimeIntegratorType, sim.mhd_timestepper)
        opts.reconstruction_type = getattr(lib.ReconstructionType, sim.reconstruction)
        opts.slope_limiter_type = getattr(lib.SlopeLimiterType, sim.limiter)
        opts.riemann_solver_type = getattr(lib.RiemannSolverType, sim.riemann)
        opts.Hall = sim.hall
        opts.Resistivity = sim.res
        opts.HyperResistivity = sim.hyper_res

    return opts
