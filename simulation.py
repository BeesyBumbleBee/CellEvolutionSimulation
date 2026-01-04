from typing import Dict
import numpy as np
import matplotlib.pyplot as plt

from chemistry import Atom
from environment import Environment, SourceType
from cell import Protocell, Genome, PossibleReactant, ReactionRule


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

    @property
    def population(self) -> int:
        return len([x for x in self.cells if x.alive])

    def setup_environment(self, sources: Dict[str, int] = None):

        if sources is None:
            sources = {
                'energy': 2,
                'co2': 3,
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
        num_reactants = int(self.rng.integers(low=2, high=4))
        num_rules = int(self.rng.integers(low=1, high=6))
        rules = [
            ReactionRule(
                num_reactants=num_reactants,
                possible_reactants=[PossibleReactant(self.rng.choice(list(Atom.basic_atoms.keys())))
                     for _ in range(num_reactants)],
                use_energy=int(self.rng.integers(0, 200)),
                priority=round(self.rng.uniform(low=0.5, high=2.0), 2)
            ) for _ in range(num_rules)
        ]

        return Genome(rules=sorted(rules, key=lambda r: r.priority))

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

    def reproduce_cell(self, cell: Protocell):
        new_genome = cell.genome.mutate()
        initial_energy = cell.reproduction_cost
        x = self.rng.choice([x for x in [cell.x-1, cell.x+1] if 0 < x < self.env.width])
        y = self.rng.choice([y for y in [cell.y-1, cell.y+1] if 0 < y < self.env.height])
        self.cells.append(
            Protocell(
                x=x,
                y=y,
                genome=new_genome,
                initial_energy=initial_energy,
            )
        )


    def step(self):
        self.env.step()
        for cell in self.cells:
            cell.step(self.env)
            if cell.can_reproduce:
                self.reproduce_cell(cell)


        self.time += 1

    def visualize(self):
        print(f" ====== t: {self.time} ====== ")
        env_fig, axes = self.env.visualize()

        for cell in [x for x in self.cells if x.alive]:
            for ax in axes:
                ax.plot(cell.x, cell.y, '.', markersize=1)

        plt.show()

        for cell in self.cells:
            print(cell.summary())
        print("\n\n")

    def run(self, steps: int = 100, visualize_steps: int = 0):
        for i in range(steps):
            self.step()
            if visualize_steps > 0 and (i % visualize_steps == 0 or i == steps - 1):
               self.visualize()

            if self.population <= 0:
                print(f"Extinction! at time step {self.time}")
                break

