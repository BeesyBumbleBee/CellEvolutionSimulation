from dataclasses import dataclass
from chemistry import Compound
from enum import StrEnum
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np


class SourceType(StrEnum):
    energy = "energy"
    n2 = "N2"
    co2 = "CO2"


@dataclass
class Source:
    x: int
    y: int
    intensity: float
    radius: int
    type: SourceType

    def get_value_at(self, target_x: int, target_y: int) -> float:
        distance = np.sqrt((self.x - target_x)**2 + (self.y - target_y)**2)
        if distance > self.radius:
            return 0.0

        falloff = 1.0 - (distance / self.radius)
        return self.intensity * (falloff ** 2)


class Environment:

    def _get_compound(self, formula: str, temperature: float) -> Compound:
        compound = Compound.from_formula(formula)
        compound.remaining_energy = round(temperature)
        return compound

    def __init__(self, width: int, height: int, ambient_temperature: float = 298.0):
        self.width: int = width
        self.height: int = height

        self.ambient_temperature: float = ambient_temperature
        self.temperature_grid = np.asarray([[ambient_temperature for w in range(width)] for h in range(height)], dtype=np.float64)

        self.grids: Dict[str, np.ndarray] = {
            source_type: np.zeros((height, width)) for source_type in SourceType
        }

        self.ambient_compounds: Dict[str, np.ndarray] = {
            'H2O': np.ones((height, width))
        }

        self.sources: List[Source] = []

        self.time_step = 0

    @property
    def ambient_temperature_grid(self) -> np.ndarray:
        return np.asarray([[self.ambient_temperature for w in range(self.width)] for h in range(self.height)], dtype=np.float64)

    def step(self):
        self._update_grid()
        self.diffuse_grids(0.9)
        self._update_temperature_grid()
        self.time_step += 1

    def get_energy_at(self, x: int, y: int) -> float:
        if 0 <= x < self.width and 0 <= y < self.height:
            return self.grids[SourceType.energy][y, x]
        return 0.0

    def extract_energy_from(self, x: int, y: int, percent: float = 1.0) -> float:
        extracted = self.get_energy_at(x,y) * percent
        self.grids[SourceType.energy][x, y] -= extracted
        return extracted

    def get_compounds_at(self, x:int, y:int) -> Dict[str, float]:
        if 0 > x or x >= self.width or 0 > y or y >= self.height:
            return {}

        compounds = {
            compound: val[x, y] for compound, val in self.grids.items() if compound != SourceType.energy
        }

        for ambient_compound in self.ambient_compounds.keys():
            compounds[ambient_compound] = self.ambient_compounds[ambient_compound][x, y]

        return compounds

    def extract_compounds_from(self, x:int, y:int, percent:float) -> Dict[str, Tuple[Compound, float]]:
        if 0 > percent:
            percent = 0.0
        if percent > 1.0:
            percent = 1.0

        compounds = self.get_compounds_at(x, y)
        extracted_compounds = {}
        for compound_formula in compounds.keys():
            if compound_formula in self.ambient_compounds.keys():
                continue
            compounds[compound_formula] *= percent
            self.grids[compound_formula][x, y] -= compounds[compound_formula]
            extracted_compounds[compound_formula] = (self._get_compound(compound_formula, self.get_temperature_at(x, y)), compounds[compound_formula])
        return extracted_compounds

    def get_temperature_at(self, x: int, y: int):
        if 0 <= x < self.width and 0 <= y < self.height:
            return self.temperature_grid[y, x]
        return 0.0

    def diffuse_grids(self, diffusion_rate: float = 0.1):
        """Simple diffusion - compounds spread to neighboring cells"""
        for source_type, grid in self.grids.items():
            new_grid = grid.copy()

            for y in range(self.height):
                for x in range(self.width):
                    if grid[y, x] > 0:
                        amount = grid[y, x] * diffusion_rate
                        neighbors = []
                        if x > 0: neighbors.append((y, x - 1))
                        if x < self.width - 1: neighbors.append((y, x + 1))
                        if y > 0: neighbors.append((y - 1, x))
                        if y < self.height - 1: neighbors.append((y + 1, x))

                        if neighbors:
                            per_neighbor = amount / len(neighbors)
                            new_grid[y, x] -= amount
                            for ny, nx in neighbors:
                                new_grid[ny, nx] += per_neighbor

            self.grids[source_type] = new_grid

    def add_source(self, x:int, y:int, intensity:float, radius: int, source_type: SourceType = SourceType.energy):
        source = Source(x, y, intensity, radius, source_type)
        self.sources.append(source)
        self._update_grid()
        self._update_temperature_grid()

    def _update_grid(self):
        for source in self.sources:
            self.grids[source.type] += np.asarray([[source.get_value_at(x, y) for x in range(self.width)] for y in range(self.height)])

    def _update_temperature_grid(self):
        self.temperature_grid = np.add(self.ambient_temperature_grid, self.grids[SourceType.energy] // 100)

    def visualize(self):
        """Create visualization of environment state"""
        fig, axes = plt.subplots(len(self.grids)+1, figsize=(4, len(self.grids)*3))

        for i, grid in enumerate(self.grids.items()):
            source_type, grid = grid
            im1 = axes[i].imshow(grid, cmap='hot', interpolation='nearest')
            axes[i].set_title(f'{source_type.title()} Distribution (t={self.time_step})')
            axes[i].set_xlabel('X Position')
            axes[i].set_ylabel('Y Position')
            #plt.colorbar(im1, ax=axes, label=f'{source_type.title()}')

            for source in [src for src in self.sources if src.type == source_type]:
                axes[i].plot(source.x, source.y, 'x', markersize=4)
                # circle = plt.Circle((source.x, source.y), source.radius,
                #                     fill=False, color='blue', linestyle='--')
                # axes[i].add_patch(circle)

        im1 = axes[-1].imshow(self.temperature_grid, cmap='hot', interpolation='nearest')
        axes[-1].set_title(f'Temperature (t={self.time_step})')
        axes[-1].set_xlabel('X Position')
        axes[-1].set_ylabel('Y Position')
        #plt.colorbar(im1, ax=axes, label=f'Temperature [K]')

        fig.tight_layout()

        return fig


if __name__ == "__main__":
    env = Environment(width=64, height=64)
    env.add_source(x=np.random.randint(20,50), y=np.random.randint(20,50), intensity=np.random.uniform(10.0, 200.0), radius=np.random.randint(1, 3), source_type=SourceType.energy)
    env.add_source(x=10, y=10, intensity=0.5, radius=1, source_type=SourceType.co2)

    for x in range(500):
        env.step()
        compounds_extracted = env.extract_compounds_from(5, 5, 1.0)
        print(compounds_extracted)
        if x % 100 == 0:
            fig = env.visualize()
            plt.show()
    fig = env.visualize()
    plt.show()

    print()
