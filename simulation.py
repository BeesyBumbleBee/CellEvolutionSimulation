from typing import Dict
import numpy as np
import matplotlib.pyplot as plt

from environment import Environment, SourceType
from cell import Protocell, Genome

class Simulation:
    def __init__(self, width: int = 32, height: int = 32):
        self.env = Environment(width, height)
        self.cells = []
        self.time = 0

        self.history = {
            'time': [],
            'population': [],
        }

        self.rng = np.random.default_rng()

    def setup_environment(self, sources: Dict[str, int] = None):

        if sources is None:
            sources = {
                'energy': 4,
                'co2': 4,
                'n2': 3,
            }

        for source_type, source_number in sources.items():
            source_type = SourceType[source_type]
            for _ in range(source_number):
                self.env.add_source(
                    x=int(self.rng.integers(low=self.env.width//10, high=(9*self.env.width)//10)),
                    y=int(self.rng.integers(low=self.env.height//10, high=(9*self.env.height)//10)),
                    intensity=self.rng.poisson(lam=(4.0 if source_type == SourceType.energy else 2.0)),
                    radius=int(self.rng.integers(low=1, high=4)),
                    source_type=source_type,
                )

    def random_genome(self):
        return Genome([])

    def populate_environment(self, num_cells: int = 10):
        for _ in range(num_cells):
            self.cells.append(
                Protocell(
                    x=int(self.rng.integers(low=0, high=self.env.width)),
                    y=int(self.rng.integers(low=0, high=self.env.height)),
                    genome=self.random_genome(),
                    initial_energy=self.rng.normal(loc=1000.0, scale=50.0),
                    )
            )


    def step(self):
        self.env.step()
        for cell in self.cells:
            cell.step(self.env)

    def visualize(self):
        env_fig, axes = self.env.visualize()

        for cell in [x for x in self.cells if x.alive]:
            for ax in axes:
                ax.plot(cell.x, cell.y, '.', markersize=1)

        plt.show()

        for cell in self.cells:
            print(cell)
        print("\n\n")

    def run(self, steps: int = 100, visualize_steps: int = 0):
        for i in range(steps):
            self.step()
            if visualize_steps > 0 and (i % visualize_steps == 0 or i == steps - 1):
               self.visualize()

