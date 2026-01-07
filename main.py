from random import randint
from os import makedirs, path

from simulation import SimulationOneShot, SimulationEpochs
from environment import Environment
from chemistry import Atom
import datetime
import numpy as np
from numpy.random import default_rng, random


def test_one_shot_sim():
    timesteps = 400
    early_stop = 600

    seed: int = np.random.default_rng().integers(low=0, high=1000000)
    ended_successfully = False

    Environment.natural_resources = [
        'energy',
        'H2O1',
        'H2',
        'O2',
        'N2',
        'Cl2',
        'C1O2',
        'P1',
        'S1',
    ]

    Atom.basic_atoms = {
        'H': {'mass': 2, 'symbol': 'H', 'electrons_in_covalence': 1, 'optimal_electrons': 2},
        'C': {'mass': 12, 'symbol': 'C', 'electrons_in_covalence': 4, 'optimal_electrons': 8},
        'N': {'mass': 14, 'symbol': 'N', 'electrons_in_covalence': 5, 'optimal_electrons': 8},
        'O': {'mass': 16, 'symbol': 'O', 'electrons_in_covalence': 6, 'optimal_electrons': 8},
        'S': {'mass': 32, 'symbol': 'S', 'electrons_in_covalence': 6, 'optimal_electrons': 8},
        'P': {'mass': 31, 'symbol': 'P', 'electrons_in_covalence': 5, 'optimal_electrons': 8},
        # 'Si': {'mass': 28, 'symbol': 'Si', 'electrons_in_covalence': 4, 'optimal_electrons': 8},
        # 'F': {'mass': 19, 'symbol': 'F', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        'Cl': {'mass': 35, 'symbol': 'Cl', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        # 'Br': {'mass': 80, 'symbol': 'Br', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        # 'I': {'mass': 127, 'symbol': 'I', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
    }


    while not ended_successfully:
        sim = SimulationOneShot()
        sim.rng = np.random.default_rng(seed=seed)
        print(f"Seed: {seed}")
        seed += 1
        sources = {
            'energy': 1,
            'S1': 2,
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
        sim.setup_environment(sources, ambient, 32, 32)
        sim.populate_environment(50)
        ended_successfully = sim.run(timesteps, timesteps // 10, early_stop_condition=early_stop)

    save_dir = f'simulation results - seed {seed-1}'
    print(f"Seed: {seed-1}")
    makedirs(save_dir, exist_ok=True)
    with open(path.join(save_dir, 'simulation_summary.txt'), 'w') as f:
        f.write(sim.summary())
    sim.graph_all(save_dir=save_dir)


def test_generations_sim():
    timesteps_per_epoch = 240
    early_stop = 300
    epochs = 15
    visualize = False

    seed: int = np.random.default_rng().integers(low=0, high=1000000)
    ended_successfully = False

    Environment.natural_resources = [
        'energy',
        'H2O1',
        'H2S1',
        'O2',
        'N1H3',
        'Cl1H1',
        'C1O2',
        'H3P1O4',
    ]

    Atom.basic_atoms = {
        'H': {'mass': 2, 'symbol': 'H', 'electrons_in_covalence': 1, 'optimal_electrons': 2},
        'C': {'mass': 12, 'symbol': 'C', 'electrons_in_covalence': 4, 'optimal_electrons': 8},
        'N': {'mass': 14, 'symbol': 'N', 'electrons_in_covalence': 5, 'optimal_electrons': 8},
        'O': {'mass': 16, 'symbol': 'O', 'electrons_in_covalence': 6, 'optimal_electrons': 8},
        'S': {'mass': 32, 'symbol': 'S', 'electrons_in_covalence': 6, 'optimal_electrons': 8},
        'P': {'mass': 31, 'symbol': 'P', 'electrons_in_covalence': 5, 'optimal_electrons': 8},
        # 'Si': {'mass': 28, 'symbol': 'Si', 'electrons_in_covalence': 4, 'optimal_electrons': 8},
        # 'F': {'mass': 19, 'symbol': 'F', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        'Cl': {'mass': 35, 'symbol': 'Cl', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        # 'Br': {'mass': 80, 'symbol': 'Br', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
        # 'I': {'mass': 127, 'symbol': 'I', 'electrons_in_covalence': 7, 'optimal_electrons': 8},
    }

    while not ended_successfully:
        sim = SimulationEpochs()
        sim.rng = np.random.default_rng(seed=seed)
        print(f"\n{datetime.datetime.now().time()} Seed: {seed}")
        seed += 1
        sources = {
            'energy': 4,
            'C1O2': 5,
            'Cl1H1': 1,
            'O2': 1,
            'N1H3': 2,
            'H2S1': 2,
            'H3P1O4': 1,
        }
        ambient = {
            'energy': 8,
            'H2O1': 5,
            'C1O2': 1,
            'O2': 1,
            'N1H3': 1,
        }
        sim.ambient = ambient
        sim.sources = sources
        sim.setup_environment(sources, ambient, 48, 48)
        sim.populate_environment(100)
        ended_successfully = sim.run(timesteps_per_epoch, epochs, visualize=visualize, early_stop_condition=early_stop)

    save_dir = f'epoch simulation results - seed {seed-1}'
    print(f"Seed: {seed - 1}")
    makedirs(save_dir, exist_ok=True)
    with open(path.join(save_dir, 'simulation_summary.txt'), 'w') as f:
        f.write(sim.summary())
    sim.graph_all(save_dir=save_dir, show=False)


if __name__ == "__main__":
    # test_one_shot_sim()
    test_generations_sim()