import sys
from pathlib import Path
from matplotlib import pyplot as plt

import pyphare.core.box as boxm
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare


atol = 1e-12

outputpath = Path("phare_outputs/comapare_diags/")
outputs = str(outputpath)


def make_fig(hier, fig_name, ilvl, collections):
    hier.plot_2d_patches(0, collections=collections).savefig(fig_name + ".png")


def plot_file(time, qty):
    return outputs + f"/{qty}_{time}.png"


def main():
    run0 = Run(sys.argv[1])
    run1 = Run(sys.argv[2])

    idx = 0

    times = run0.all_times()["B"]
    print("times", times)
    for time in times:
        print("\ntime :", time)

        # eqr = hierarchy_compare(
        #     run0.GetB(time, all_primal=False),
        #     run1.GetB(time, all_primal=False),
        #     atol=atol,
        # )
        # print("\n\t B", eqr)

        # e0 = run0.GetE(time, all_primal=False)
        # e1 = run1.GetE(time, all_primal=False)
        # eqr = hierarchy_compare(e0, e1, atol=atol)
        # print("\n\t E", eqr)

        # if not eqr:
        #     ediff = VectorField(e0) - VectorField(e1)
        #     for c in ["x", "y", "z"]:
        #         ediff.plot(
        #             filename=f"E{c}_diff_{time}.png",
        #             qty=f"{c}",
        #             plot_patches=True,
        #             vmin=0,
        #             vmax=atol,
        #         )

        N0 = run0.GetNi(time)
        N1 = run1.GetNi(time)
        eqr = hierarchy_compare(N0, N1, atol=atol)
        print(f"\n\n Ni", eqr)
        # if not eqr:
        Ndiff = ScalarField(N0) - ScalarField(N1)
        Ndiff.plot(
            filename=plot_file(time, "Ni"),
            plot_patches=True,
            vmin=0,
            vmax=atol,
        )

        pops = run0.all_pops()
        for pop in pops:
            N0 = run0.GetN(time, pop_name=pop)
            N1 = run1.GetN(time, pop_name=pop)
            eqr = hierarchy_compare(N0, N1, atol=atol)
            print(f"\n\n rho {pop}", eqr)
            # if not eqr:
            Ndiff = ScalarField(N0) - ScalarField(N1)
            Ndiff.plot(
                filename=plot_file(time, "N"),
                plot_patches=True,
                vmin=0,
                vmax=atol,
            )

            # eqr = hierarchy_compare(
            #     run0.GetParticles(time, pop), run1.GetParticles(time, pop), atol=atol
            # )
            # print(f"\n\n domain particles {pop}", eqr)

            # f0 = run0.GetFlux(time, pop_name=pop)
            # f1 = run1.GetFlux(time, pop_name=pop)
            # eqr = hierarchy_compare(f0, f1, atol=atol)
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
    outputpath.mkdir(parents=True, exist_ok=True)
    main()
