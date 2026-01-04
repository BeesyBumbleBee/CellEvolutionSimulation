from simulation import Simulation

def test_sim():
    sim = Simulation(8, 8)
    sources = {
        'energy': 1,
        'co2': 2,
        'n2': 1,
    }
    sim.setup_environment(sources)
    sim.populate_environment(1)

    sim.run(400, 50)



if __name__ == "__main__":
    test_sim()