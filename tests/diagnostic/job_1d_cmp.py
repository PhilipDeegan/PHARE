import traceback

print("!cmp!")
try:
    import pyphare.pharein as ph
    from pyphare.pharein import ElectronModel
    from tests.simulator import basicSimulatorArgs, makeBasicModel
    from tests.diagnostic import dump_all_diags

    out = "phare_outputs/diags_1d/cmp"
    simInput = {
        "refinement_boxes": {},
        "diag_options": {
            "format": "phareh5",
            "options": {"dir": out, "mode": "overwrite"},
        },
    }

    ph.Simulation(**basicSimulatorArgs(ndim=1, interp=1, **simInput))
    model = makeBasicModel()
    ElectronModel(closure="isothermal", Te=0.12)
    dump_all_diags(model.populations)

except Exception:
    print('Exception caught in "Simulator.setup()": {}'.format(sys.exc_info()))
    print(traceback.format_exc())
    raise ValueError("Error in Simulator.setup(), see previous error")
