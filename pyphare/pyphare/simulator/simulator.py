#
#
#


import atexit
import time as timem
import numpy as np
import pyphare.pharein as ph


life_cycles = {}


@atexit.register
def simulator_shutdown():
    import gc
    from ._simulator import obj

    if obj is not None:  # needs to be killed before MPI
        obj.reset()
    life_cycles.clear()
    gc.collect()  # force cause we need to be sure it's unloaded


# see src/core/data/particles/particle_array_def.hpp
def make_cpp_simulator(
    cpp_lib, dim, interp, nbrRefinedPart, hier, layout=1, allocator=0
):
    make_sim = f"make_simulator_{dim}_{interp}_{nbrRefinedPart}_{layout}_{allocator}"
    print("make_sim", make_sim)
    return getattr(cpp_lib, make_sim)(hier)


def startMPI():
    from pyphare.cpp import cpp_lib

    if "samrai" not in life_cycles:
        life_cycles["samrai"] = cpp_lib().SamraiLifeCycle()


class Simulator:
    def __init__(self, simulation, auto_dump=True, **kwargs):
        assert isinstance(simulation, ph.Simulation)  # pylint: disable=no-member
        self.simulation = simulation
        self.cpp_hier = None  # HERE
        self.cpp_sim = None  # BE
        self.cpp_dw = None  # DRAGONS, i.e. use weakrefs if you have to ref these.
        self.post_advance = kwargs.get("post_advance", None)

        self.print_eol = "\n"
        if kwargs.get("print_one_line", True):
            self.print_eol = "\r"
        self.print_eol = kwargs.get("print_eol", self.print_eol)
        self.log_to_file = kwargs.get("log_to_file", True)

        self.auto_dump = auto_dump
        import pyphare.simulator._simulator as _simulator

        _simulator.obj = self
        from pyphare.cpp import cpp_lib

        self.cpp_lib = cpp_lib()

    def __del__(self):
        self.reset()

    def setup(self, **kwargs):
        # mostly to detach C++ class construction/dict parsing from C++ Simulator::init
        try:
            import pyphare.cpp.validate as validate_cpp

            startMPI()

            if all([not self.simulation.dry_run, self.simulation.write_reports]):
                # not necessary during testing
                validate_cpp.log_runtime_config()
            validate_cpp.check_build_config_is_runtime_compatible()

            if self.log_to_file:
                self._log_to_file()
            ph.populateDict()
            self.cpp_hier = self.cpp_lib.make_hierarchy()
            self.cpp_sim = make_cpp_simulator(
                self.cpp_lib,
                self.simulation.ndim,
                self.simulation.interp_order,
                self.simulation.refined_particle_nbr,
                self.cpp_hier,
                **kwargs,
            )
            return self
        except:
            import sys

            print(
                'Exception caught in "Simulator.setup()": {}'.format(sys.exc_info()[0])
            )
            raise ValueError("Error in Simulator.setup(), see previous error")

    def initialize(self):
        try:
            if self.cpp_hier is None:
                self.setup()

            if self.simulation.dry_run:
                return self

            self.cpp_sim.initialize()
            self._auto_dump()  # first dump might be before first advance
            return self
        except:
            import sys

            print(
                'Exception caught in "Simulator.initialize()": {}'.format(
                    sys.exc_info()[0]
                )
            )
            raise ValueError("Error in Simulator.initialize(), see previous error")
        print("initialize ")

    def _throw(self, e):
        import sys

        if self.cpp_lib.mpi_rank() == 0:
            print(e)
        sys.exit(1)

    def advance(self, dt=None):
        self._check_init()
        if self.simulation.dry_run:
            return self
        if dt is None:
            dt = self.timeStep()

        try:
            self.cpp_sim.advance(dt)
        except (RuntimeError, TypeError, NameError, ValueError) as e:
            self._throw(f"Exception caught in simulator.py::advance: \n{e}")
        except KeyboardInterrupt as e:
            self._throw(f"KeyboardInterrupt in simulator.py::advance: \n{e}")

        if self._auto_dump() and self.post_advance != None:
            self.post_advance(self.cpp_sim.currentTime())
        return self

    def times(self):
        return np.arange(
            self.cpp_sim.startTime(),
            self.cpp_sim.endTime() + self.timeStep(),
            self.timeStep(),
        )

    def run(self):
        self._check_init()
        if self.simulation.dry_run:
            return self
        perf = []
        end_time = self.cpp_sim.endTime()
        t = self.cpp_sim.currentTime()

        while t < end_time:
            tick = timem.time()
            self.advance()
            tock = timem.time()
            ticktock = tock - tick
            perf.append(ticktock)
            t = self.cpp_sim.currentTime()
            if self.cpp_lib.mpi_rank() == 0:
                out = f"t = {t:8.5f}  -  {ticktock:6.5f}sec  - total {np.sum(perf):7.4}sec"
                print(out, end=self.print_eol)

        print("mean advance time = {}".format(np.mean(perf)))
        print("total advance time = {}".format(np.sum(perf)))

        return self.reset()

    def _auto_dump(self):
        return self.auto_dump and self.dump()

    def dump(self, *args):
        assert len(args) == 0 or len(args) == 2

        if len(args) == 0:
            return self.cpp_sim.dump(
                timestamp=self.currentTime(), timestep=self.timeStep()
            )
        return self.cpp_sim.dump(timestamp=args[0], timestep=args[1])

    def data_wrangler(self):
        self._check_init()
        if self.cpp_dw is None:
            from pyphare.data.wrangler import DataWrangler

            self.cpp_dw = DataWrangler(self)
        return self.cpp_dw

    def reset(self):
        if self.cpp_sim is not None:
            ph.clearDict()
        if self.cpp_dw is not None:
            self.cpp_dw.kill()
        self.cpp_dw = None
        self.cpp_sim = None
        self.cpp_hier = None
        if "samrai" in life_cycles:
            type(life_cycles["samrai"]).reset()
        return self

    def timeStep(self):
        self._check_init()
        return self.cpp_sim.timeStep()

    def currentTime(self):
        self._check_init()
        return self.cpp_sim.currentTime()

    def domain_box(self):
        self._check_init()
        return self.cpp_sim.domain_box()

    def cell_width(self):
        self._check_init()
        return self.cpp_sim.cell_width()

    def interp_order(self):
        self._check_init()
        return self.cpp_sim.interp_order  # constexpr static value

    def _check_init(self):
        if self.cpp_sim is None:
            self.initialize()

    def _log_to_file(self):
        """
        send C++ std::cout logs to files with env var PHARE_LOG
        Support keys:
            RANK_FILES - logfile per rank
            DATETIME_FILES - logfile with starting datetime timestamp per rank
            NONE - no logging files, display to cout
        """
        import os

        if "PHARE_LOG" not in os.environ:
            os.environ["PHARE_LOG"] = "RANK_FILES"

        if os.environ["PHARE_LOG"] != "NONE" and self.cpp_lib.mpi_rank() == 0:
            from pathlib import Path

            Path(".log").mkdir(exist_ok=True)
