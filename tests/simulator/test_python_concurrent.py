"""
  This script exists to minimize testing time by running all simulation/tests
    concurrently without needing to wait for any particular file or set of tests
"""

import os
import unittest
import multiprocessing

from tests.simulator.test_validation import SimulatorValidation

from tests.simulator.initialize.test_fields_init_1d import InitializationTest as InitField1d
from tests.simulator.initialize.test_particles_init_1d import InitializationTest as InitParticles1d
from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
from tests.simulator.advance.test_particles_advance_1d import AdvanceTest as AdvanceParticles1d

from tests.simulator.initialize.test_fields_init_2d import InitializationTest as InitField2d
from tests.simulator.initialize.test_particles_init_2d import InitializationTest as InitParticles2d
from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
from tests.simulator.advance.test_particles_advance_2d import AdvanceTest as AdvanceParticles2d

from tests.simulator.initialize.test_fields_init_3d import InitializationTest as InitField3d
from tests.simulator.initialize.test_particles_init_3d import InitializationTest as InitParticles3d
from tests.simulator.advance.test_fields_advance_3d import AdvanceTest as AdvanceField3d
from tests.simulator.advance.test_particles_advance_3d import AdvanceTest as AdvanceParticles3d


N_CORES = int(os.environ["N_CORES"]) if "N_CORES" in os.environ else multiprocessing.cpu_count()
DO_DIM = int(os.environ["DO_DIM"]) if "DO_DIM" in os.environ else 0
MPI_RUN = int(os.environ["MPI_RUN"]) if "MPI_RUN" in os.environ else 1
PRINT   = int(os.environ["PRINT"]) if "PRINT" in os.environ else 0
ALL_DIMS = DO_DIM == 0

def test_cmd(clazz, test_id):
    return f"mpirun -n {MPI_RUN} python3 -m {clazz.__module__} {clazz.__name__}.{test_id}"

if __name__ == "__main__":

    test_classes_to_run = []

    if ALL_DIMS:
        test_classes_to_run += [SimulatorValidation]

    if ALL_DIMS or DO_DIM == 1:
        test_classes_to_run += [
            InitField1d,
            InitParticles1d,
            AdvanceField1d,
            AdvanceParticles1d,
        ]
    if ALL_DIMS or DO_DIM == 2:
        test_classes_to_run += [
            InitField2d,
            InitParticles2d,
            AdvanceField2d,
            AdvanceParticles2d,
        ]
    if ALL_DIMS or DO_DIM == 3:
        test_classes_to_run += [
            InitField3d,
            InitParticles3d,
            AdvanceField3d,
            AdvanceParticles3d
        ]

    tests = []
    loader = unittest.TestLoader()
    for test_class in test_classes_to_run:
        for suite in loader.loadTestsFromTestCase(test_class):
            tests += [test_cmd(type(suite), suite._testMethodName)]

    if PRINT:
        # handy to get a list of all individual test cases as command line strings
        for test in tests:
            print(test)
    else:
        from tools.python3 import run_mp
        try:
            run_mp(tests, N_CORES, check=True)
        except:
            import sys
            sys.exit(1)
