import sys

import pyphare.core.box as boxm
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from matplotlib import pyplot as plt
from tests.simulator import diff_boxes


def make_fig(hier, fig_name, ilvl, collections):
    hier.plot_2d_patches(0, collections=collections).savefig(fig_name + ".png")


def fail_fig(hier, eqr):
    for res, ref, cmp in eqr.failed:
        error_boxes = diff_boxes(ref[:], cmp[:], ref.box, atol=1e-15)

        print("ref.box", ref.box)

        extra_collections = [
            {
                "boxes": error_boxes,
                "facecolor": "black",
            }
        ]
        make_fig(hier, ref.name, 0, extra_collections)


def main():
    run0 = Run(sys.argv[1])
    run1 = Run(sys.argv[2])
    atol = 1e-15
    idx = 0
    #
    for time in run0.all_times()["B"]:
        print("\ntime :", time)

        pops = run0.all_pops()
        print("pops", pops)
        for pop in pops:
            eqr = hierarchy_compare(
                run0.GetParticles(time, pop), run1.GetParticles(time, pop), atol=atol
            )
            print(f"\n\n domain particles {pop}", eqr)

            eqr = hierarchy_compare(
                run0.GetParticles(time, pop, kind="levelGhost"),
                run1.GetParticles(time, pop, kind="levelGhost"),
                atol=atol,
            )
            print(f"\n\n level ghost particles {pop}", eqr)

        eqr = hierarchy_compare(
            run0.GetB(time, all_primal=False),
            run1.GetB(time, all_primal=False),
            atol=atol,
        )
        print("\n\t B", eqr)

        e0 = run0.GetE(time, all_primal=False)
        e1 = run1.GetE(time, all_primal=False)
        eqr = hierarchy_compare(e0, e1, atol=atol)
        print("\n\t E", eqr)
        # fail_fig(e0, eqr)
        # (e0 - e1).plot(qty="x", filename="ediff.png", plot_patches=True)

        # pops = run0.all_pops()
        # print("pops", pops)
        # for pop in pops:
        #     eqr = hierarchy_compare(
        #         run0.GetParticles(time, pop), run1.GetParticles(time, pop), atol=atol
        #     )
        #     print(f"\n\n domain particles {pop}", eqr)

        # eqr = hierarchy_compare(
        #     run0.GetN(time, pop_name=pop), run1.GetN(time, pop_name=pop), atol=atol
        # )
        # print(f"\n\n rho {pop}", eqr)

        # f0 = run0.GetFlux(time, pop_name=pop)
        # f1 = run1.GetFlux(time, pop_name=pop)
        # eqr = hierarchy_compare(f0, f1, atol=atol)
        # fail_fig(f0, eqr)
        # print(f"\n\n F {pop}", eqr)

        # eqr = hierarchy_compare(run0.GetVi(time), run1.GetVi(time), atol=atol)
        # print(f"\n\n GetVi ", eqr)
        # ni0 = run0.GetNi(time)

        # fig, ax = plt.subplots()
        # ni0.plot(
        #     ax=ax,
        #     qty="value",
        #     # filename="ni_master",
        #     plot_patches=True,
        #     vmin=0.25,
        #     vmax=1.5,
        #     marker="+",
        # )
        # ni1 = run1.GetNi(time).plot(
        #     ax=ax,
        #     qty="value",
        #     # filename="ni_lol",
        #     plot_patches=True,
        #     vmin=0.25,
        #     vmax=1.5,
        #     marker="o",
        #     ls="-",
        # )
        # fig.savefig("lol.png")
        # eqr = hierarchy_compare(run0.GetNi(time), run1.GetNi(time), atol=atol)
        # print(f"\n\n GetNi ", eqr)


if __name__ == "__main__":
    main()
