
import sys
import unittest
import numpy as np
from ddt import ddt, data, unpack
from tests.diagnostic import all_timestamps
import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, parse_cli_args

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

from pyphare.pharein.simulation import supported_dimensions

permutations = [
    [dim, res] for dim in supported_dimensions() for res in ["low", "high"]
]

configs = {
    1 : {
        "modes0" : [4, 8, 16, 32, 64, 128, 256, 512],
        "modes1" : [4, 8, 16, 32, 64, 128, 256, 512],
        "b_amplitudes0" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        "b_amplitudes1" : [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005],
        "low" :{
            "time_step_nbr": 300000,
            "cells"        : 2000,
            "dl"           : 0.2,
        },
        "high":{
            "time_step_nbr": 800000,
            "cells"        : 4000,
            "dl"           : 0.1,
        },
        "base":{ # this is overridden by resolution specific keys if they exist
            "smallest_patch_size":20,
            "largest_patch_size":50,
            "final_time":200.,
            "boundary_types":"periodic",
            "diag_options":{"format": "phareh5",
                          "options": {"dir": "setOfModes1d",
                                      "mode":"overwrite"}}
        },
    },
    2 : {
        "modes0" : [4, 8, 16, 32, 64],
        "modes1" : [4, 8, 16, 32, 64],
        "b_amplitudes0" : [0.002, 0.002, 0.002, 0.002, 0.002],
        "b_amplitudes1" : [0.005, 0.005, 0.005, 0.005, 0.005],
        "low" :{
            "time_step_nbr": 150000,
            "cells"        : (400, 20),
            "dl"           : (0.2, 0.2),
        },
        "high":{
            "time_step_nbr": 800000,
            "cells"        : (200, 10),
            "dl"           : (0.1, 0.1),
        },
        "base":{ # this is overridden by resolution specific keys if they exist
            "smallest_patch_size":10,
            "largest_patch_size":10,
            "final_time":200.,
            "boundary_types":("periodic", "periodic"),
            "diag_options":{"format": "phareh5",
                          "options": {"dir": "setOfModes2d",
                                      "mode":"overwrite"}}
        },
    },
}

def get_theta():
    from pyphare.pharein.global_vars import sim
    L = sim.cells
    return np.arctan2(L[1], L[0])


def setup(dim, resolution, modes, b_amplitudes, polarization, seed=12345):
    assert(len(modes) == len(b_amplitudes))

    ph.global_vars.sim = None

    config = configs[dim]["base"].copy()
    config.update(configs[dim][resolution])
    sim = ph.Simulation(**config)

    # list of wave_numbers for the given box
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()
    wave_numbers = [2*np.pi*m/L[0] for m in modes]

    if dim == 2:
        theta = get_theta()
        wave_num_x = [k*np.cos(theta) for k in wave_numbers]
        wave_num_y = [k*np.sin(theta) for k in wave_numbers]


    all_funcs = {
        "base" :{ # this is overridden by dimension specific keys if they exist
            "density": lambda *xyz: 1,
            "vthx": lambda *xyz: .01,
            "vthy": lambda *xyz: .01,
            "vthz": lambda *xyz: .01,
            "vbulkx": lambda *xyz: 0,
            "vbulky": lambda *xyz: 0,
            "vbulkz": lambda *xyz: 0,
        },
        1: {
            "bx": lambda *xyz: 1,
            "by": lambda x: np.sum([b*np.cos(k*x) for (k, b) in zip(wave_numbers, b_amplitudes)]),
            "bz": lambda x: np.sum([b*np.sin(k*x)*polarization for (k, b) in zip(wave_numbers, b_amplitudes)]),
        },
        2: { # could be moved to a function see, https://github.com/PHAREHUB/PHARE/blob/master/tests/simulator/__init__.py#L93 # rm comment
            "bx": lambda x, y: np.cos(theta) - np.sum([b*np.cos(kx*x+ky*y)*np.sin(theta) for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
            "by": lambda x, y: np.sin(theta) + np.sum([b*np.cos(kx*x+ky*y)*np.cos(theta) for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
            "bz": lambda x, y: np.sum([b*np.sin(kx*x+ky*y)*polarization for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
        },
    }

    funcs = all_funcs["base"].copy()
    funcs.update(all_funcs[dim])

    vvv = {
        "vbulkx": funcs["vbulkx"], "vbulky": funcs["vbulky"], "vbulkz": funcs["vbulkz"],
        "vthx": funcs["vthx"], "vthy": funcs["vthy"], "vthz": funcs["vthz"]
    }
    ph.MaxwellianFluidModel(
        bx=funcs["bx"], by=funcs["by"], bz=funcs["bz"],
        main={"charge": 1, "density": funcs["density"], **vvv, "init":{"seed": seed}}
    )
    ph.ElectronModel(closure="isothermal", Te=0.)

    timestamps = all_timestamps(sim)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )

    return sim, wave_numbers, b_amplitudes


def post_sim0(dim, **kwargs):
    if cpp.mpi_rank() == 0 and dim == 2:
        theta = get_theta()
    if cpp.mpi_rank() == 0:
        # do things
        pass
    cpp.mpi_barrier() # KEEP THIS!


def post_sim1(dim, **kwargs):
    if cpp.mpi_rank() == 0 and dim == 2:
        theta = get_theta()
    if cpp.mpi_rank() == 0:
        # do things
        pass
    cpp.mpi_barrier() # KEEP THIS!



def dispersion(dim, resolution):
    seed = cpp.mpi_rank()+1
    config = configs[dim][resolution]

    modes = configs[dim]["modes0"]
    b_amplitudes = configs[dim]["b_amplitudes0"]
    polarization = -1
    sim, wave_nums, b1 = setup(dim, resolution, modes, b_amplitudes, polarization, seed)
    Simulator(sim).run().reset()
    post_sim0(sim)

    modes = configs[dim]["modes1"]
    b_amplitudes = configs[dim]["b_amplitudes1"]
    polarization = 1
    sim, wave_nums, b1 = setup(dim, resolution, modes, b_amplitudes, polarization, seed)
    Simulator(sim).run().reset()
    post_sim1(sim)



class DispersionTest(SimulatorTest):
    # ddt mangles method names so direct lookup isn't easy
    def test_dispersion(self, dim, resolution):
        dispersion(dim, resolution)

@ddt
class DDTDispersionTest(SimulatorTest):
    @data(*permutations)
    @unpack
    def test_dispersion(self, dim, resolution):
        dispersion(dim, resolution)


if __name__=="__main__":

    if len(sys.argv) == 3:
        dim, resolution = parse_cli_args()
        if resolution not in ['low', 'high']:
            raise ValueError('arg should be "low" or "high"')

        DispersionTest().setUp().test_dispersion(int(dim), resolution).tearDown()

    elif len(sys.argv) == 1:

        loader = unittest.TestLoader()
        suites = []
        for suite in loader.loadTestsFromTestCase(DDTDispersionTest):
            suites += [suite]
        tests = unittest.TestSuite(suites)
        unittest.TextTestRunner(verbosity=2).run(tests)

    else:
        print('example usage: $script $dim $resolution=[low/high]')

