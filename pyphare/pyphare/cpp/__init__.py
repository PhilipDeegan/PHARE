#
#
#


import json

# continue to use override if set
_cpp_lib_override = None


def cpp_lib(override=None):
    import importlib

    def _try(mod, fallthrough=True):
        try:
            return importlib.import_module(f"pybindlibs.{mod}")
        except ImportError as e:
            print(e)
            if not fallthrough:
                raise e

    global _cpp_lib_override

    if override is not None:
        _cpp_lib_override = override

    if _cpp_lib_override is not None:
        return _try(_cpp_lib_override)

    if not __debug__:
        return importlib.import_module("pybindlibs.cpp")
    try:
        return importlib.import_module("pybindlibs.cpp_dbg")
    except ImportError:
        return importlib.import_module("pybindlibs.cpp")


def cpp_etc_lib():
    import importlib

    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(dim, interp, n_particles):
    return getattr(cpp_lib(), f"Splitter_{dim}_{interp}_{n_particles}")


def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def split_pyarrays_fn(dim, interp, n_particles):
    return getattr(cpp_lib(), f"split_pyarray_particles_{dim}_{interp}_{n_particles}")


def supported_particle_layouts():
    # see: src/core/data/particles/particle_array_def.hpp

    return [1]  # , 3
