
import numpy as np

def splitter_type(dim, interp, n_particles):
    from pybindlibs import cpp
    return getattr(cpp,
        "Splitter_" + str(dim) + "_" + str(interp)+ "_" + str(n_particles),
    )

def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()



def cpp_float_id(dtype):

    if dtype in [np.float32, np.dtype("float32"), "float32"]:
        return "s"

    if dtype in [np.float64, np.dtype("float64"), "float64"]:
        return "d"

    raise ValueError("Unsupported dtype")



def split_pyarrays_fn(dim, interp, n_particles, dtype=np.float64):
    from pybindlibs import cpp
    return getattr(cpp,
        "split_pyarray_particles_" + str(dim) + "_" + str(interp)+ "_" + str(n_particles) + "_" + cpp_float_id(dtype),
    )

