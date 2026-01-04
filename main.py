from simulation import Simulation
from numpy.random import default_rng
import matplotlib.pyplot as plt


def test_sim():
    end_time = 0
    timesteps = 200
    seed: int = 379

    while end_time < timesteps:
        sim = Simulation(8, 8)
        sim.rng = default_rng(seed=seed)
        seed += 1
        sources = {
            'energy': 0,
            'n1': 1,
            'h2': 1,
            'c1': 2,
        }
        sim.setup_environment(sources)
        sim.populate_environment(5)

        sim.run(timesteps, timesteps // 10)
        end_time = sim.time

    labels = list(sim.history.keys())
    values = list(sim.history.values())
    fig = plt.scatter(x=values[0], y=values[1])
    plt.xlabel('Time')
    plt.ylabel('Population')
    print(f"Seed: {seed-1}")
    plt.show()
    pass


if __name__ == "__main__":
    test_sim()