from environment import Environment, SourceType
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt

class Simulation:
    def __init__(self, width: int = 32, height: int = 32):
        self.env = Environment(width, height)
        self.cells = []
        self.time = 0

        self.history = {
            'time': [],
            'population': [],
        }

        self.rng = np.random

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
                    x=self.rng.randint(round(self.env.width * 0.1), round(self.env.width * 0.9)),
                    y=self.rng.randint(round(self.env.height * 0.1), round(self.env.height * 0.9)),
                    intensity=self.rng.poisson(lam=(100.0 if source_type == SourceType.energy else 2.0)),
                    radius=self.rng.randint(1, 4),
                    source_type=source_type,
                )


    def step(self):
        self.env.step()

    def visualize(self):
        fig = self.env.visualize()

        plt.show()

    def run(self, steps: int = 100, visualize_steps: int = 0):
        for i in range(steps):
            self.step()
            if visualize_steps > 0 and (i % visualize_steps == 0 or i == steps - 1):
                self.visualize()

