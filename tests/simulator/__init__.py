

from pybindlibs import cpp
import pyphare.pharein as ph, numpy as np
from pyphare.pharein import ElectronModel

def basicSimulatorArgs(dim:int, interp:int, **kwargs):

    cells = [65 for i in range(dim)]
    if "cells" in kwargs:
        cells = kwargs["cells"]
    if not isinstance(cells, list):
        cells = [cells]
    dl = [1.0 / v for v in cells]
    b0 = [[10 for i in range(dim)], [50 for i in range(dim)]]
    boundary = ["periodic" for i in range(dim)]
    args = {
        "interp_order": interp,
        "smallest_patch_size": 5,
        "largest_patch_size": 64,
        "time_step_nbr": 1000,
        "final_time": 1.0,
        "boundary_types": boundary,
        "cells": cells,
        "dl": dl,
        "max_nbr_levels": 2,
        "refinement_boxes": {"L0": {"B0": b0}},
        "refined_particle_nbr": dim * 2,
        "diag_options": {},
    }
    for k, v in kwargs.items():
        if k in args:
            args[k] = v
    return args


def defaultPopulationSettings():
    background_particles = 0.1  # avoids 0 density
    xmax = ph.globals.sim.simulation_domain()[0]
    pi_over_xmax = np.pi / xmax
    pop_fn = {
        1:{
        "density": lambda x: 1.0 / np.cosh((x - xmax * 0.5)) ** 2 + background_particles,
        "vbulkx": lambda x: np.sin(1 * pi_over_xmax * x),
        "vbulky": lambda x: np.sin(1 * pi_over_xmax * x),
        "vbulkz": lambda x: np.sin(1 * pi_over_xmax * x),
        "vthx": lambda x: 1,
        "vthy": lambda x: 1,
        "vthz": lambda x: 1,
      },2 :{
        "density": lambda x, y: 1.0 / np.cosh((x - xmax * 0.5)) ** 2 + background_particles,
        "vbulkx": lambda x, y: np.sin(1 * pi_over_xmax * x),
        "vbulky": lambda x, y: np.sin(1 * pi_over_xmax * x),
        "vbulkz": lambda x, y: np.sin(1 * pi_over_xmax * x),
        "vthx": lambda x, y: 1,
        "vthy": lambda x, y: 1,
        "vthz": lambda x, y: 1,
      }
    }
    dim = len(ph.globals.sim.cells)
    return {
        "charge": 1,
        "density": pop_fn[dim]["density"],
        "vbulkx": pop_fn[dim]["vbulkx"],
        "vbulky": pop_fn[dim]["vbulky"],
        "vbulkz": pop_fn[dim]["vbulkz"],
        "vthx": pop_fn[dim]["vthx"],
        "vthy": pop_fn[dim]["vthy"],
        "vthz": pop_fn[dim]["vthz"],
    }


def makeBasicModel(extra_pops={}):
    xmax = ph.globals.sim.simulation_domain()[0]
    pi_over_xmax = np.pi / xmax
    EM_fn = {
        1:{
        "bx":lambda x: np.cos(2 * pi_over_xmax * x),
        "by":lambda x: np.sin(1 * pi_over_xmax * x),
        "bz":lambda x: np.cos(2 * pi_over_xmax * x),
      },2 :{
        "bx":lambda x, y: np.cos(2 * pi_over_xmax * x),
        "by":lambda x, y: np.sin(1 * pi_over_xmax * x),
        "bz":lambda x, y: np.cos(2 * pi_over_xmax * x),
      }
    }
    pops = {
        "protons": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
        "alpha": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 1000,
            "init": {"seed": 13337},
        },
    }
    pops.update(extra_pops)
    dim = len(ph.globals.sim.cells)
    return ph.MaxwellianFluidModel(
        bx=EM_fn[dim]["bx"],
        by=EM_fn[dim]["by"],
        bz=EM_fn[dim]["bz"],
        **pops
    )



def create_simulator(dim, interp, **input):
    cpp.reset()
    ph.globals.sim = None
    ph.Simulation(**basicSimulatorArgs(dim, interp, **input))
    extra_pops = {}
    if "populations" in input:
        for pop, vals in input["populations"].items():
            extra_pops[pop] = defaultPopulationSettings()
            extra_pops[pop].update(vals)

    model = makeBasicModel(extra_pops)
    if "diags_fn" in input:
        input["diags_fn"](model)

    ElectronModel(closure="isothermal",Te = 0.12)

    ph.populateDict()
    import sys
    try:
        hier = cpp.make_hierarchy()
        sim = cpp.make_simulator(hier)
        sim.initialize()
        return [cpp.make_diagnostic_manager(sim, hier), sim, hier]
    except:
        e = sys.exc_info()[0]
        print('Exception caught in "create_simulator": {}'.format(e))
