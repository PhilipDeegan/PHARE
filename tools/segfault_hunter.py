

import tools.bench.real.bench_harris as harris

if __name__=="__main__":

    import random
    from pyphare.simulator.simulator import Simulator

    while True:
        rando = random.randint(0, 1e10)
        harris.seed = rando = random.randint(0, 1e10)
        print("seed", harris.seed)
        harris.config()
        Simulator(ph.global_vars.sim).run()
