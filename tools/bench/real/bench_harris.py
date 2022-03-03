#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics,FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import MetaDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein import global_vars as gv
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from pyphare.cpp import cpp_lib
cpp = cpp_lib()
startMPI()


diag_outputs="phare_outputs/harris"

def config():
    L=0.5
    Simulation(
        smallest_patch_size=4,
        time_step=0.01,
        final_time=40.,
        #boundary_types="periodic",
        cells=(100,100),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels = 3,
        nesting_buffer=2,
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": diag_outputs,
                                  "mode":"overwrite"}}
    )

    def density(x, y):
        from pyphare.pharein.global_vars import sim
        Ly = sim.simulation_domain()[1]
        return 0.4 + 1./np.cosh((y-Ly*0.3)/L)**2 + 1./np.cosh((y-Ly*0.7)/L)**2


    def S(y, y0, l):
        return 0.5*(1. + np.tanh((y-y0)/l))


    def by(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - 0.5 * Lx)
        y1 = (y - 0.3 * Ly)
        y2 = (y - 0.7 * Ly)

        dBy1 = 2*dB*x0 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBy2 = -2*dB*x0 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        return dBy1 + dBy2


    def bx(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - 0.5 * Lx)
        y1 = (y - 0.3 * Ly)
        y2 = (y - 0.7 * Ly)

        dBx1 = -2*dB*y1 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBx2 =  2*dB*y2 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(y,Ly*0.3,L) -S(y, Ly*0.7,L)) + dBx1 + dBx2


    def bz(x, y):
        return 0.


    def b2(x, y):
        return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2


    def T(x, y):
        K = 0.7
        temp = 1./density(x, y)*(K - b2(x, y)*0.5)
        assert np.all(temp >0)
        return temp

    def vx(x, y):
        return 0.


    def vy(x, y):
        return 0.


    def vz(x, y):
        return 0.


    def vthx(x, y):
        return np.sqrt(T(x, y))


    def vthy(x, y):
        return np.sqrt(T(x, y))


    def vthz(x, y):
        return np.sqrt(T(x, y))


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz,
        "nbr_part_per_cell":100
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,  **vvv}#, "init":{"seed":cpp.mpi_rank()+int(sys.argv[1])}}
    )

    ElectronModel(closure="isothermal", Te=0.0)



    sim = ph.global_vars.sim
    dt =   1.*sim.time_step
    nt = sim.final_time/dt+1
    timestamps = dt * np.arange(nt)
    print(timestamps)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )

    MetaDiagnostics(
        quantity="tags",
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
    )
  # for popname in ("protons",):
  #     for name in ["domain", "levelGhost", "patchGhost"]:
  #         ParticleDiagnostics(quantity=name,
  #                             compute_timestamps=timestamps,
  #                             write_timestamps=timestamps,
  #                             population_name=popname)



def get_time(path, time, datahier = None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from
    datahier = hierarchy_from(h5_filename=path+"/EM_E.h5", time=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path+"/EM_B.h5", time=time, hier=datahier)
    return datahier

def post_advance(time):
    if cpp.mpi_rank() == 0:
        print(f"running tests at time {time}")

        for lvl_idx,level in get_time(diag_outputs, time).levels().items():
            for patch_i, ref_patch in enumerate(level.patches):
                for ref_pdname, ref_pd in ref_patch.patch_datas.items():
                    for cmp_patch in level.patches[patch_i + 1:]:
                        print("ref_pdname", ref_pdname, lvl_idx)
                        cmp_pd = cmp_patch.patch_datas[ref_pdname]
                        if cmp_pd.box * ref_pd.box is not None:
                            raise ValueError(f"Domains overlap {cmp_pd.box}, {ref_pd.box}")

        print(f"tests passed")
    cpp.mpi_barrier()

def main():

    config()
    Simulator(gv.sim, post_advance=post_advance).run()


if __name__=="__main__":
    main()






# import numpy as np
# from pyphare.cpp import cpp_lib # must be first
# cpp_lib("pybindlibs.cpp_sim_2_1_4")
# import pyphare.pharein as ph
# from pyphare.core.box import nDBox

# seed = 133333333337
# cells, dl = 400, .2
# patch_sizes = [20,20]
# time_step_nbr=100
# time_step=0.001
# diag_outputs="tools/bench/real/harris/outputs"

# def density(x, y):
#     L = ph.global_vars.sim.simulation_domain()[1]
#     return 0.2 + 1./np.cosh((y-L*0.3)/0.5)**2 + 1./np.cosh((y-L*0.7)/0.5)**2

# def by(x, y):
#     sim = ph.global_vars.sim
#     Lx = sim.simulation_domain()[0]
#     Ly = sim.simulation_domain()[1]
#     w1, w2 = 0.2, 1.0
#     x0 = (x - 0.5 * Lx)
#     y1 = (y - 0.3 * Ly)
#     y2 = (y - 0.7 * Ly)
#     w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
#     w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
#     w5 = 2.0*w1/w2
#     return (w5 * x0 * w3) + ( -w5 * x0 * w4)

# def S(y, y0, l): return 0.5*(1. + np.tanh((y-y0)/l))
# def bx(x, y):
#     sim = ph.global_vars.sim
#     Lx = sim.simulation_domain()[0]
#     Ly = sim.simulation_domain()[1]
#     w1, w2 = 0.2, 1.0
#     x0 = (x - 0.5 * Lx)
#     y1 = (y - 0.3 * Ly)
#     y2 = (y - 0.7 * Ly)
#     w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
#     w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
#     w5 = 2.0*w1/w2
#     v1, v2 = -1, 1.
#     return v1 + (v2-v1)*(S(y,Ly*0.3,0.5) -S(y, Ly*0.7, 0.5)) + (-w5*y1*w3) + (+w5*y2*w4)

# def bz(x, y): return 0.
# def b2(x, y): return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2
# def T(x, y):  return 1./density(x, y)*(1 - b2(x, y)*0.5)
# def vxyz(x, y): return 0.
# def vthxyz(x, y): return np.sqrt(T(x, y))

# def config():
#     ph.Simulation(# strict=True,
#         smallest_patch_size=patch_sizes[0], largest_patch_size=patch_sizes[1],
#         time_step_nbr=time_step_nbr, time_step=time_step,
#         cells=[cells] * 2, dl=[dl] * 2,
#         resistivity=0.001, hyper_resistivity=0.001,
#         diag_options={"format": "phareh5", "options": {"dir": diag_outputs, "mode":"overwrite"}},
#         refinement="tagging", max_nbr_levels=2, nesting_buffer=2,
#     )
#     ph.MaxwellianFluidModel( bx=bx, by=by, bz=bz,
#         protons={"charge": 1, "density": density, "init":{"seed": seed},
#           **{ "nbr_part_per_cell":100,
#             "vbulkx": vxyz, "vbulky": vxyz, "vbulkz": vxyz,
#             "vthx": vthxyz, "vthy": vthxyz, "vthz": vthxyz,
#           }
#         },
#     )
#     ph.ElectronModel(closure="isothermal", Te=0.0)

#     from tests.diagnostic import all_timestamps
#     timestamps = all_timestamps(ph.global_vars.sim)
#     timestamps = np.asarray([timestamps[0], timestamps[-1]])
#     for quantity in ["E", "B"]:
#         ph.ElectromagDiagnostics(
#             quantity=quantity,
#             write_timestamps=timestamps,
#             compute_timestamps=timestamps,
#         )

# if ph.PHARE_EXE or __name__=="__main__":
#     config()

# if __name__=="__main__":
#     from pyphare.simulator.simulator import Simulator
#     Simulator(ph.global_vars.sim).run()
