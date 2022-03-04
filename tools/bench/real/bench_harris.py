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
from pyphare.pharesee.run import Run
cpp = cpp_lib()
startMPI()

time_step=0.01
diag_outputs="phare_outputs/harris"
dl=.4
L0_to_L2_ratio = 4

def config():
    L=0.5
    Simulation(
        smallest_patch_size=4,
        time_step=time_step,
        time_step_nbr=100,
        #final_time=40.,
        #boundary_types="periodic",
        cells=(100,100),
        dl=(dl, dl),
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


def plot(time, ref_pd, cmp_pd):
    r = Run(diag_outputs)
    Ni = r.GetNi(time)

    import numpy as np
    b0 = ref_pd.box
    b1 = cmp_pd.box
    lower_x = np.min([b0.lower[0], b1.lower[0]])
    upper_x = np.max([b0.upper[0], b1.upper[0]])
    lower_y = np.min([b0.lower[1], b1.lower[0]])
    upper_y = np.max([b0.upper[1], b1.upper[1]])
    x_lim = [lower_x / L0_to_L2_ratio * dl - 5, upper_x / L0_to_L2_ratio * dl + 5]
    y_lim = [lower_y / L0_to_L2_ratio * dl - 5, upper_y / L0_to_L2_ratio * dl + 5]

    fig, ax = Ni.plot(qty="rho", vmin=-0.1, vmax=3.6, plot_patches=True,
        patchcolors=['k','r','k'], ls=["-", "-", "-"],levels=(2,),
       lw=[1.5,0.5,0.8], xlim=x_lim, ylim=y_lim)
    fig.savefig(f"ni_{time}.png")


def get_time(path, time, datahier = None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from
    datahier = hierarchy_from(h5_filename=path+"/EM_B.h5", time=time, hier=datahier)
    return datahier

def post_advance(time):
    if cpp.mpi_rank() == 0:
        print(f"running tests at time {time}")
        for lvl_idx,level in get_time(diag_outputs, time).levels().items():
            for patch_i, ref_patch in enumerate(level.patches):
                for ref_pdname, ref_pd in ref_patch.patch_datas.items():
                    for cmp_patch in level.patches[patch_i + 1:]:
                        cmp_pd = cmp_patch.patch_datas[ref_pdname]
                        if cmp_pd.box * ref_pd.box is not None:
                            plot(time, ref_pd, cmp_pd)
                            raise ValueError(f"Domains overlap {cmp_pd.box}, {ref_pd.box}")
    cpp.mpi_barrier()

def main():
    config()
    Simulator(gv.sim, post_advance=post_advance).run()

if __name__=="__main__":
    main()

    # for i in range(55, 66):
    #     time = time_step * i
    #     post_advance(time)

