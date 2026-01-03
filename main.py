from simulation import Simulation

def test_sim():
    sim = Simulation(32, 32)
    sources = {
        'energy': 1,
        'co2': 2,
        'n2': 2,
    }
    sim.setup_environment(sources)

    sim.run(100, 10)



if __name__ == "__main__":
    test_sim()