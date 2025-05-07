#
#
#


import json

# continue to use override if set
_cpp_lib_override = None
_cpp_lib = None


def cpp_lib(override=None):
    global _cpp_lib
    if _cpp_lib:
        return _cpp_lib

    import importlib

    global _cpp_lib_override
    if override is not None:
        _cpp_lib_override = override
    if _cpp_lib_override is not None:
        _cpp_lib = importlib.import_module(_cpp_lib_override)
        return _cpp_lib

    if not __debug__:
        _cpp_lib = importlib.import_module("pybindlibs.cpp")
    try:
        _cpp_lib = importlib.import_module("pybindlibs.cpp_dbg")
    except ImportError:
        _cpp_lib = importlib.import_module("pybindlibs.cpp")
    return _cpp_lib


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
