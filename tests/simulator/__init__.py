import os
import unittest
from datetime import datetime

import numpy as np
import pyphare.pharein as ph
from pyphare.core.phare_utilities import np_array_ify


def parse_cli_args(pop_from_sys=True):
    import sys

    r = sys.argv[1:].copy()
    if pop_from_sys:  # args can interfere with other things
        sys.argv = [sys.argv[0]]
    return r


# Block accidental dictionary key rewrites
class NoOverwriteDict(dict):
    def __init__(self, dict):
        for k, v in dict.items():
            self[k] = v

    def __setitem__(self, k, v):
        if k in self.keys():
            raise ValueError("Key is already present")
        else:
            return super(NoOverwriteDict, self).__setitem__(k, v)


def basicSimulatorArgs(ndim: int, interp: int, **kwargs):
    from pyphare.pharein.simulation import valid_refined_particle_nbr
    from pyphare.pharein.simulation import check_patch_size
    from pyphare.core.phare_utilities import np_array_ify

    cells = np_array_ify(kwargs.get("cells", 20), ndim)

    _, smallest_patch_size = check_patch_size(ndim, interp_order=interp, cells=cells)
    dl = 1.0 / cells
    b0 = [[3] * ndim, [8] * ndim]

    args = {
        "interp_order": interp,
        "smallest_patch_size": smallest_patch_size,
        "largest_patch_size": [20] * ndim,
        "time_step_nbr": 1000,
        "final_time": 1.0,
        "boundary_types": ["periodic"] * ndim,
        "cells": cells,
        "dl": dl,
        "refinement_boxes": {"L0": {"B0": b0}},
        "refined_particle_nbr": valid_refined_particle_nbr[ndim][interp][0],
        "diag_options": {},
        "nesting_buffer": 0,
        "strict": True,
    }
    for k, v in kwargs.items():
        if k in args:
            args[k] = v
    args["cells"] = np_array_ify(args["cells"], ndim)
    args["dl"] = np_array_ify(args["dl"], ndim)
    return args


def meshify(*xyz):
    if all([isinstance(v, np.ndarray) for v in xyz]):
        return xyz
    return np.meshgrid(*xyz, indexing="ij")


def pi_over_max_domain():
    return [np.pi / max_domain for max_domain in ph.global_vars.sim.simulation_domain()]


def fn_periodic(sim, *xyz):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()
    _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
    return np.asarray([_(i) for i, v in enumerate(xyz)]).prod(axis=0)


def density_1d_periodic(sim, x):
    xmax = sim.simulation_domain()[0]
    background_particles = 0.3  # avoids 0 density
    return 1.0 / np.cosh((x - xmax * 0.5) ** 2 + background_particles)


def density_2d_periodic(sim, x, y):
    xmax, ymax = sim.simulation_domain()
    background_particles = 0.3  # avoids 0 density
    xx, yy = meshify(x, y)
    return (
        np.exp(-((xx - 0.5 * xmax) ** 2)) * np.exp(-((yy - ymax / 2.0) ** 2))
        + background_particles
    )


# def density_3d_periodic(sim, x, y, z):
#     xmax, ymax, zmax = sim.simulation_domain()
#     background_particles = 0.3  # avoids 0 density
#     xx, yy, zz = meshify(x, y, z)
#     r = np.exp(-(xx-0.5*xmax)**2)*np.exp(-(yy-ymax/2.)**2)*np.exp(-(zz-zmax/2.)**2) + background_particles
#     return r


def defaultPopulationSettings(sim, density_fn, vbulk_fn):
    return {
        "charge": 1,
        "density": lambda *xyz: density_fn(sim, *xyz),
        "vbulkx": lambda *xyz: vbulk_fn(sim, *xyz) + 0.01,
        "vbulky": lambda *xyz: vbulk_fn(sim, *xyz) + 0.02,
        "vbulkz": lambda *xyz: vbulk_fn(sim, *xyz) + 0.03,
        "vthx": lambda *xyz: 1,
        "vthy": lambda *xyz: 1,
        "vthz": lambda *xyz: 1,
    }


def makeBasicModel(extra_pops={}, ppc=100):
    sim = ph.global_vars.sim
    _density_fn_periodic = globals()["density_" + str(sim.ndim) + "d_periodic"]

    pops = {
        "protons": {
            **defaultPopulationSettings(sim, _density_fn_periodic, fn_periodic),
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
        "alpha": {
            **defaultPopulationSettings(sim, _density_fn_periodic, fn_periodic),
            "nbr_part_per_cell": ppc,
            "init": {"seed": 13337},
        },
    }
    pops.update(extra_pops)
    return ph.MaxwellianFluidModel(
        bx=lambda *xyz: fn_periodic(sim, *xyz) + 0.04,
        by=lambda *xyz: fn_periodic(sim, *xyz) + 0.05,
        bz=lambda *xyz: fn_periodic(sim, *xyz) + 0.06,
        **pops,
    )


def populate_simulation(dim, interp, **input):
    ph.global_vars.sim = None
    simulation = ph.Simulation(**basicSimulatorArgs(dim, interp, **input))
    extra_pops = {}
    if "populations" in input:
        for pop, vals in input["populations"].items():
            extra_pops[pop] = defaultPopulationSettings()
            extra_pops[pop].update(vals)

    model = makeBasicModel(extra_pops)
    if "diags_fn" in input:
        input["diags_fn"](model)

    ph.ElectronModel(closure="isothermal", Te=0.12)

    return simulation


def diff_boxes(slice1, slice2, box, atol=None):
    from pyphare.core.box import Box

    if atol is not None:
        ignore = np.isclose(slice1, slice2, atol=atol, rtol=0)

        def _diff(slice0):
            slice0[ignore] = 0  # set values which are within atol range to 0
            return slice0

        diff = np.abs(_diff(slice1.copy()) - _diff(slice2.copy()))
    else:
        diff = np.abs(slice1 - slice2)

    boxes = []
    if box.ndim == 1:
        x1 = np.where(diff != 0)
        for x in zip(x1):
            x = x + box.lower[0]
            boxes += [Box([x], [x])]
    elif box.ndim == 2:
        x1, y1 = np.where(diff != 0)
        for x, y in zip(x1, y1):
            x = x + box.lower[0]
            y = y + box.lower[1]
            boxes += [Box([x, y], [x, y])]
    elif box.ndim == 3:
        x1, y1, z1 = np.where(diff != 0)
        for x, y, z in zip(x1, y1, z1):
            x = x + box.lower[0]
            y = y + box.lower[1]
            z = z + box.lower[2]
            boxes += [Box([x, y, z], [x, y, z])]
    return boxes


#
#


def caliper_func_times_json(data_dir, mpi_rank=0):
    return f"{os.path.join(data_dir, f'func_times.{mpi_rank}.json')}"


def caliper_recorder_cali(data_dir, mpi_rank=0):
    return f"{os.path.join(data_dir, f'recorder.{mpi_rank}.cali')}"


CALIPER_MODES = [
    # "callpath:event:recorder:trace",
    "report,event,trace,timestamp,recorder",  # light
    "alloc,aggregate,cpuinfo,memusage,debug,env,event,loop_monitor,region_monitor,textlog,io,pthread,sysalloc,recorder,report,timestamp,statistics,spot,trace,validator,mpi,mpireport,mpiflush",  # heavy
]


def activate_caliper(data_dir, mode_idx=0):
    from pyphare.cpp import cpp_lib

    rank = cpp_lib().mpi_rank()
    env = os.environ

    # env["CALI_SERVICES_ENABLE"] = "event,trace,timer,report"
    # env["CALI_REPORT_CONFIG"] = "format json"
    # env["CALI_REPORT_FILENAME"] = "trace.json"

    # env[
    #     "CALI_CONFIG"
    # ] = "hatchet-region-profile,topdown-counters.all,output.format=json"

    # # env["CALI_CONFIG_PROFILE"] = "callstack-trace"
    # env["CALI_SERVICES_ENABLE"] = CALIPER_MODES[mode_idx]

    # env["CALI_CONFIG"] = "hatchet-region-profile"

    # # env["CALI_CALLPATH_USE_NAME"] = "true"

    # env["CALI_REPORT_FILENAME"] = caliper_func_times_json(data_dir, rank)
    # env[
    #     "CALI_REPORT_CONFIG"
    # ] = "SELECT function,time.duration ORDER BY time.duration FORMAT json"
    # env["CALI_RECORDER_FILENAME"] = caliper_recorder_cali(data_dir, rank)

    # print("os.environ", os.environ)


class SimulatorTest(unittest.TestCase):
    test_kwargs = ["rethrow"]

    def tearDown(self):
        self.clean_up_diags_dirs()

    def setUp(self):
        from pyphare.simulator.simulator import startMPI

        startMPI()

    def datetime_now(self):
        return datetime.now()

    def datetime_diff(self, then):
        return (datetime.now() - then).total_seconds()

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    def pop(kwargs):
        for key in SimulatorTest.test_kwargs:
            if key in kwargs:
                kwargs.pop(key)
        return kwargs

    old_failureException = unittest.TestCase.failureException

    @property  # intercept test failure to not delete diags in case
    def failureException(self):
        self.success = False
        return self.old_failureException

    def register_diag_dir_for_cleanup(self, diag_dir):
        self.diag_dirs += [diag_dir]

    def __init__(self, *args, **kwargs):
        super(SimulatorTest, self).__init__(*args, **SimulatorTest.pop(kwargs.copy()))
        self.rethrow_ = True
        for key in SimulatorTest.test_kwargs:
            if key in kwargs:
                super().__setattr__(f"{key}_", kwargs[key])
        self.diag_dirs = []  # cleanup after tests
        self.success = True

    def run(self, result=None):
        self._outcome = result
        super().run(result)

    def unique_diag_dir_for_test_case(self, base_path, ndim, interp, post_path=""):
        from pyphare.cpp import cpp_lib

        cpp = cpp_lib()
        return f"{base_path}/{self._testMethodName}/{cpp.mpi_size()}/{ndim}/{interp}/{post_path}"

    def clean_up_diags_dirs(self):
        from pyphare.cpp import cpp_lib

        cpp_lib().mpi_barrier()
        if cpp_lib().mpi_rank() == 0 and self.success:
            import os
            import shutil

            for diag_dir in self.diag_dirs:
                if os.path.exists(diag_dir):
                    shutil.rmtree(diag_dir)
        cpp_lib().mpi_barrier()
