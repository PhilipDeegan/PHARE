
objects  = {}
diag_info = {}
sim = None

def fn_type():
    import os
    init = "PHARE_INIT_FN"
    if init in os.environ:
        return os.environ[init]
    return "vector"

func_type = fn_type()
