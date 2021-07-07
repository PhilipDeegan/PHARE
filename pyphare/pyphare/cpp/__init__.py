
import numpy as np

def cpp_lib():
    import importlib
    if not __debug__:
        return importlib.import_module("pybindlibs.cpp")
    try:
        return importlib.import_module("pybindlibs.cpp_dbg")
    except ImportError as err:
        return importlib.import_module("pybindlibs.cpp")


def splitter_type(dim, interp, n_particles):
    return getattr(cpp_lib(), f"Splitter_{dim}_{interp}_{n_particles}")

def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def cpp_float_id(dtype):

    if dtype in [np.float32, np.dtype("float32"), "float32"]:
        return "s"

    if dtype in [np.float64, np.dtype("float64"), "float64"]:
        return "d"

    raise ValueError("Unsupported dtype")

def split_pyarrays_fn(dim, interp, n_particles, dtype=np.float64):
    return getattr(cpp_lib(), f"split_pyarray_particles_{dim}_{interp}_{n_particles}_{cpp_float_id(dtype)}")


