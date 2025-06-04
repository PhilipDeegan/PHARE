import sys

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from matplotlib import pyplot as plt

run0 = Run(sys.argv[1])
run1 = Run(sys.argv[2])
atol = 1e-15
idx = 0
#
for time in run0.all_times()["B"]:
    print("\ntime :", time)

    eqr = hierarchy_compare(
        run0.GetB(time, all_primal=False), run1.GetB(time, all_primal=False), atol=atol
    )
    print("\n\t B", eqr)

    eqr = hierarchy_compare(
        run0.GetE(time, all_primal=False), run1.GetE(time, all_primal=False), atol=atol
    )
    print("\n\t E", eqr)

    pops = run0.all_pops()
    for pop in pops:
        eqr = hierarchy_compare(
            run0.GetParticles(time, pop_name=pop),
            run1.GetParticles(time, pop_name=pop),
            atol=atol,
        )
        print(f"\n\t domain particles {pop}", eqr)

        eqr = hierarchy_compare(
            run0.GetN(time, pop_name=pop), run1.GetN(time, pop_name=pop), atol=atol
        )
        print(f"\n\t rho {pop}", eqr)

        eqr = hierarchy_compare(
            run0.GetFlux(time, pop_name=pop),
            run1.GetFlux(time, pop_name=pop),
            atol=atol,
        )
        print(f"\n\t F {pop}", eqr)

        eqr = hierarchy_compare(
            run0.GetParticles(time, pop), run1.GetParticles(time, pop), atol=atol
        )
        print(f"\n\t Particles {pop}", eqr)

    eqr = hierarchy_compare(run0.GetVi(time), run1.GetVi(time), atol=atol)
    print(f"\n\t GetVi ", eqr)
    ni0 = run0.GetNi(time)

    fig, ax = plt.subplots()
    ni0.plot(
        ax=ax,
        qty="value",
        # filename="ni_master",
        plot_patches=True,
        vmin=0.25,
        vmax=1.5,
        marker="+",
    )
    ni1 = run1.GetNi(time).plot(
        ax=ax,
        qty="value",
        # filename="ni_lol",
        plot_patches=True,
        vmin=0.25,
        vmax=1.5,
        marker="o",
        ls="-",
    )
    fig.savefig("lol.png")
    eqr = hierarchy_compare(run0.GetNi(time), run1.GetNi(time), atol=atol)
    print(f"\n\t GetNi ", eqr)
