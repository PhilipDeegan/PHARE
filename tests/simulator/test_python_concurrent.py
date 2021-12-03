"""
  This script exists to minimize testing time by running all simulation/tests
    concurrently without needing to wait for any particular file or set of tests
"""

import os
import unittest
import multiprocessing
import concurrent.futures
from tools.python3 import run, RunTimer

from tests.simulator.test_validation import SimulatorValidation
from tests.simulator.test_diagnostics import DiagnosticsTest  # mpirun -n 1/2/3/4
from tests.simulator.initialize.test_fields_init_1d import InitializationTest as InitField1d
from tests.simulator.initialize.test_particles_init_1d import InitializationTest as InitParticles1d
from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
from tests.simulator.advance.test_particles_advance_1d import AdvanceTest as AdvanceParticles1d
from tests.simulator.initialize.test_fields_init_2d import InitializationTest as InitField2d
from tests.simulator.initialize.test_particles_init_2d import InitializationTest as InitParticles2d
from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
from tests.simulator.advance.test_particles_advance_2d import AdvanceTest as AdvanceParticles2d

N_CORES = int(os.environ["N_CORES"]) if "N_CORES" in os.environ else multiprocessing.cpu_count()
MPI_RUN = os.environ["MPI_RUN"] if "MPI_RUN" in os.environ else 1
PRINT   = int(os.environ["PRINT"]) if "PRINT" in os.environ else 0

def test_cmd(clazz, test_id, mpi_run):
    return f"mpirun -n {mpi_run} python3 -m {clazz.__module__} {clazz.__name__}.{test_id}"

test_classes_to_run = [
  SimulatorValidation, DiagnosticsTest,
  InitField1d,       InitParticles1d,
  AdvanceField1d,    AdvanceParticles1d,
  InitField2d,       InitParticles2d,
  AdvanceField2d,    AdvanceParticles2d
]

class TestBatch:
    def __init__(self, tests, mpi_run = 1):
        self.tests = tests
        self.mpi_run = mpi_run

def load_test_cases_in(classes, mpi_run = 1):
    tests, loader = [], unittest.TestLoader()
    for test_class in classes:
        for suite in loader.loadTestsFromTestCase(test_class):
            tests += [test_cmd(type(suite), suite._testMethodName, mpi_run)]
    return TestBatch(tests, mpi_run)

def build_batches():
    batches = []
    if MPI_RUN=="cmake":
        batches += [load_test_cases_in(test_classes_to_run, 1)]
        batches += [load_test_cases_in(test_classes_to_run, 2)]
        batches += [load_test_cases_in([DiagnosticsTest], 3)]
        batches += [load_test_cases_in([DiagnosticsTest], 4)]
    else:
        batches += [load_test_cases_in(test_classes_to_run), int(MPI_RUN)]
    return batches

def print_tests(batches):
    for batch in batches:
        for test in batch.tests:
            print(test)

class CallableTest:
    def __init__(self, bi, ti, cmd):
        self.bi = bi
        self.ti = ti
        self.cmd = cmd
        self.run = None

    def __call__(self, **kwargs):
        #print("Launching", self.cmd)
        self.run = RunTimer(self.cmd, shell=True, capture_output=True, check=True, print_cmd=False)
        return self


def pick_tests(batches, cc, tests_run, all_from=None):
    tests = []
    #print("1 pick_tests", len(tests), cc.cores_avail)
    if all_from == None:
        for bi, batch in enumerate(batches):
            for ti, test in enumerate(batch.tests):
                if batch.mpi_run <= cc.cores_avail:
                    if ti not in tests_run[bi]:
                        tests += [CallableTest(bi, ti, test)]
                        cc.cores_avail -= batch.mpi_run
                        tests_run[bi].add(ti)
    else:
        bi = all_from
        for ti, test in enumerate(batches[bi].tests):
            tests += [CallableTest(bi, ti, test)]
            cc.cores_avail -= batches[bi].mpi_run
            #print("- pick_tests", batches[bi].mpi_run)
            tests_run[bi].add(ti)
    #print("2 pick_tests", len(tests), cc.cores_avail)
    return tests, cc, tests_run

class CoreCount:
  def __init__(self, cores_avail):
      self.cores_avail = cores_avail

def submit(executor, batches, cc, tests_run, tests = []):
    if len(tests) == 0:
        tests, cc, tests_run = pick_tests(batches, cc, tests_run)

    if len(tests) == 0: return cc, tests_run
    jobs = [executor.submit(test) for test in tests]
    #print("submitted qSize=%i"%executor._work_queue.qsize(), cc.cores_avail)
    for future in concurrent.futures.as_completed(jobs):
        try:
            proc = future.result()
            if future.exception() is not None:
                raise future.exception()
            print(proc.run.cmd, f"finished in {proc.run.t:.2f} seconds")
            cc.cores_avail += batches[proc.bi].mpi_run
            #print("qSize=%i"%executor._work_queue.qsize(), cc.cores_avail)
            if executor._work_queue.qsize() == 0 and cc.cores_avail > 1:
                cc, tests_run = submit(executor, batches, cc, tests_run)
        except Exception as exc:
            executor.shutdown(wait=False, cancel_futures=True)
            raise exc
    return cc, tests_run

def run_tests(batches):
    cc = CoreCount(N_CORES)
    assert cc.cores_avail > max([batch.mpi_run for batch in batches])

    tests_run = [set() for batch in batches]
    with concurrent.futures.ThreadPoolExecutor(max_workers=N_CORES) as executor:
        tests, cc, tests_run = pick_tests(batches, cc, tests_run, all_from=0)
        #print("start", cc.cores_avail)
        submit(executor, batches, cc, tests_run, tests)


if __name__ == "__main__":
    batches = build_batches()
    if PRINT:
        print_tests(batches)
    else:
        run_tests(batches)
