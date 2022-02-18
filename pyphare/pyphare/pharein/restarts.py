
from ..core import phare_utilities
from . import global_vars

# ------------------------------------------------------------------------------



def restarts_checker(func):
    def wrapper(restarts_object, **kwargs):

        mandatory_keywords = ['write_timestamps']

        # check if some mandatory keywords are not missing
        missing_mandatory_kwds = phare_utilities.check_mandatory_keywords(mandatory_keywords, **kwargs)
        if len(missing_mandatory_kwds) > 0:
            raise RuntimeError("Error: missing mandatory parameters : " + ', '.join(missing_mandatory_kwds))

        accepted_keywords = ['path', 'flush_every']
        accepted_keywords += mandatory_keywords

        # check that all passed keywords are in the accepted keyword list
        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise RuntimeError("Error: invalid arguments - " + " ".join(wrong_kwds))

        try:
            # just take mandatory arguments from the dict
            # since if we arrived here we are sure they are there

            kwargs['path'] = kwargs.get("path", './')

            return func(restarts_object, **kwargs)

        except ValueError as msg:
            print(msg)

    return wrapper


import numpy as np
# ------------------------------------------------------------------------------
def validate_timestamps(clazz, **kwargs):
    sim = global_vars.sim

    for key in ["write_timestamps"]:
        timestamps = kwargs[key]

        if np.any(timestamps < sim.init_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({sim.init_time}))")
        if np.any(timestamps > sim.final_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))")
        if not np.all(np.diff(timestamps) >= 0):
            raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")
        if not np.all(np.abs(timestamps / sim.time_step - np.rint(timestamps/sim.time_step) < 1e-9)):
            raise RuntimeError(f"Error: {clazz}.{key} is inconsistent with simulation.time_step")





# ------------------------------------------------------------------------------

def try_cpp_dep_vers():
    try:
        from pyphare.cpp import cpp_etc_lib
        return cpp_etc_lib().phare_deps()
    except ImportError:
        return {}



class Restarts(object):

    h5_flush_never = 0
    cpp_dep_vers = try_cpp_dep_vers()

    @restarts_checker
    def __init__(self, **kwargs):

        if global_vars.sim is None:
            raise RuntimeError("A simulation must be created before adding restarts")

        # self.name = name
        self.path = kwargs['path']

        validate_timestamps(self.__class__.__name__, **kwargs)
        self.write_timestamps = kwargs['write_timestamps']

        global_vars.sim.add_restarts(self)


# ------------------------------------------------------------------------------


