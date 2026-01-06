from random import randint
from os import makedirs, path

from simulation import Simulation
import numpy as np
from numpy.random import default_rng, random
import matplotlib.pyplot as plt


def test_sim():
    timesteps = 5000
    early_stop = 10000

    seed: int = np.random.default_rng().integers(low=0, high=1000000)
    ended_successfully = False


    while not ended_successfully:
        sim = Simulation(64, 64)
        sim.rng = np.random.default_rng(seed=seed)
        print(f"Seed: {seed}")
        seed += 1
        sources = {
            'energy': 1,
            'S1': 2,
            # 'Cl2': 2,
            # 'I2': 1,
            # 'F2': 1,
            # 'Br2': 1,
            'C1O2': 3,
            'P1': 1,
        }
        ambient = {
            'H2O1': 5,
            'H2': 2,
            'O2': 1,
            'N2': 1,
            'Cl2': 1,
        }
        sim.setup_environment(sources, ambient)
        sim.populate_environment(50)
        ended_successfully = sim.run(timesteps, timesteps // 10, early_stop_condition=early_stop)

    save_dir = f'simulation results - seed {seed-1}'
    print(f"Seed: {seed-1}")
    makedirs(save_dir, exist_ok=True)
    with open(path.join(save_dir, 'simulation_summary.txt'), 'w') as f:
        f.write(sim.summary())
    sim.graph_all(save_dir=save_dir)




if __name__ == "__main__":
    test_sim()