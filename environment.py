from dataclasses import dataclass
from chemistry import Compound
from enum import StrEnum
from typing import Dict, List, Tuple
from math import floor, ceil, sqrt

import matplotlib.pyplot as plt
import numpy as np



@dataclass
class Source:
    x: int
    y: int
    intensity: float
    radius: int
    type: str
    mask: np.ndarray

    def __repr__(self):
        return f"Source of {self.type:<12s} at (x={self.x:>3d}, y={self.y:>3d}) with intensity {self.intensity:>4.2f} and radius {self.radius:>2d}"

    def get_value_at(self, target_x: int, target_y: int) -> float:
        distance = np.sqrt((self.x - target_x)**2 + (self.y - target_y)**2)
        if distance > self.radius:
            return 0.0

        falloff = 1.0 - (distance / self.radius)
        return self.intensity * (falloff ** 2)


class Environment:
    natural_resources = [
        'energy',
        'H2O1',
        'H2',
        'O2',
        'N2',
        'Cl2',
        'C1O2',
        'P1',
        # 'N1H3',
        # 'C1H4',
        # 'C1O2',
        # 'C1',
        # 'H2S1',
        # 'O2',
        # 'N2',
        # 'H2',
        # 'Cl2',
        # 'I2',
        # 'F2',
        # 'Br2',
        'S1',
        # 'P1',
    ]

    def __init__(self, width: int, height: int, ambient_temperature: float = 323.0):
        self.width: int = width
        self.height: int = height

        self.ambient_temperature: float = ambient_temperature
        self.temperature_grid = np.asarray([[ambient_temperature for w in range(width)] for h in range(height)], dtype=np.float64)

        self.grids: Dict[str, np.ndarray] = {
            resource_type: np.zeros((height, width)) for resource_type in Environment.natural_resources
        }

        self.ambient_compounds: Dict[str, int] = {
            resource_type: 0 for resource_type in Environment.natural_resources
        }

        self.compounds: Dict[str, Compound] = {
            compound_formula: Compound.from_formula(compound_formula) for compound_formula in Environment.natural_resources if compound_formula != 'energy'
        }

        self.sources: List[Source] = []

        self.time_step = 0

    def __repr__(self):
        return f"Environment of size {self.width}x{self.height}"

    def summary(self) -> str:
        out_str = f"{self}\n"
        out_str += "Resources spawning inside:\n"
        for resource in Environment.natural_resources:
            out_str += f"\t{resource}\n"
        out_str += "\nSources:\n"
        for source in self.sources:
            out_str += f"\t{source}\n"
        out_str += "\nAmbient resources:\n"
        for resource, val in self.ambient_compounds.items():
            if val == 0:
                continue
            out_str += f'\t{resource:12s} = {val:<4.2f}\n'
        return out_str

    @property
    def ambient_temperature_grid(self) -> np.ndarray:
        return np.asarray([[self.ambient_temperature for w in range(self.width)] for h in range(self.height)], dtype=np.float64)

    def set_ambient_resource(self, resource: str, amount: int):
        assert resource in Environment.natural_resources

        self.ambient_compounds[resource] = amount

    def _get_compound(self, formula: str, temperature: float) -> Compound:
        from copy import deepcopy
        compound = deepcopy(self.compounds[formula])
        compound.remaining_energy = round(temperature)
        return compound

    def step(self):
        self._update_grid()
        self.diffuse_grids(0.1)
        self._update_temperature_grid()
        self.time_step += 1

    def get_energy_at(self, x: int, y: int) -> float:
        if 0 <= x < self.width and 0 <= y < self.height:
            return float(self.grids['energy'][y, x])
        return 0.0

    def extract_energy_from(self, x: int, y: int, percent: float = 1.0) -> float:
        energy_at = self.get_energy_at(x,y)
        extracted = energy_at * percent
        if energy_at <= extracted or extracted < 0:
            self.grids['energy'][y, x] = 0
            return energy_at
        self.grids['energy'][y, x] -= extracted
        return extracted

    def get_compounds_at(self, x:int, y:int) -> Dict[str, int]:
        if 0 > x or x >= self.width or 0 > y or y >= self.height:
            return {}

        compounds = {
            compound_formula: floor(float(val[y, x])) for compound_formula, val in self.grids.items() if compound_formula != 'energy'
        }

        return compounds

    def extract_compounds_from(self, x:int, y:int, percent:float, max_compounds: Dict[str, int]) -> List[Compound]:
        if 0 > percent:
            percent = 0.0
        if percent > 1.0:
            percent = 1.0

        compounds = {
            compound_formula: (max(0, min(val - max_compounds[compound_formula], val)) if compound_formula in max_compounds.keys() else val)
            for compound_formula, val in self.get_compounds_at(x, y).items()
        }
        extracted_compounds: List[Compound] = []
        for compound_formula in compounds.keys():
            compounds[compound_formula] *= percent
            self.grids[compound_formula][y, x] -= floor(compounds[compound_formula])
            extracted_compounds.extend([
                self._get_compound(compound_formula, self.get_temperature_at(x, y)) for _ in range(floor(compounds[compound_formula]))
            ])
        return extracted_compounds

    def get_temperature_at(self, x: int, y: int):
        if 0 <= x < self.width and 0 <= y < self.height:
            return self.temperature_grid[y, x]
        return 0.0

    def diffuse_grids(self, diffusion_rate: float = 0.1):
        """Simple diffusion - compounds spread to neighboring cells"""

        for source_type, grid in self.grids.items():
            new_grid = grid.copy()

            def can_spread(x, y, val) -> bool:
                return grid[y, x] < val

            for y in range(self.height):
                for x in range(self.width):
                    if grid[y, x] > 0:
                        val = grid[y, x]
                        amount = grid[y, x] * diffusion_rate
                        neighbors = []
                        if x > 0 and can_spread(x-1, y, val): neighbors.append((y, x - 1))
                        if x < self.width - 1 and can_spread(x+1, y, val): neighbors.append((y, x + 1))
                        if y > 0 and can_spread(x, y-1, val): neighbors.append((y - 1, x))
                        if y < self.height - 1 and can_spread(x, y+1, val): neighbors.append((y + 1, x))

                        if neighbors:
                            per_neighbor = amount / len(neighbors)
                            new_grid[y, x] -= amount
                            for ny, nx in neighbors:
                                new_grid[ny, nx] += per_neighbor
            self.grids[source_type] = new_grid

    def _create_source_mask(self, x:int, y:int, intensity:float ,radius:int) -> np.ndarray:
        y_coords = np.arange(self.height)
        x_coords = np.arange(self.width)
        x_grid, y_grid = np.meshgrid(x_coords, y_coords)

        distances = np.sqrt((x_grid - x) ** 2 + (y_grid - y) ** 2)
        disc = (distances <= radius).astype(int)
        return disc * intensity

    def add_source(self, x:int, y:int, intensity:float, radius: int, resource_type: str = 'energy'):
        assert resource_type in Environment.natural_resources
        if intensity < 0.0:
            intensity = -intensity
        source = Source(x, y, intensity, radius, resource_type, self._create_source_mask(x,y,intensity,radius))
        self.sources.append(source)
        self._update_grid()
        self._update_temperature_grid()

    def _update_grid(self):
        for resource, val in self.ambient_compounds.items():
            self.grids[resource] += np.asarray(
                [[max(0, val - self.grids[resource][y, x])for x in range(self.width)] for y in range(self.height)])

        for source in self.sources:
            self.grids[source.type] += source.mask

    def _update_temperature_grid(self):
        self.temperature_grid = np.add(self.ambient_temperature_grid, self.grids['energy'] // 10)

    def visualize(self):
        """Create visualization of environment state"""
        x_len = floor(sqrt(len(self.grids) + 1))
        y_len = ceil(sqrt(len(self.grids) + 1))
        fig, axes = plt.subplots(x_len, y_len, figsize=(6*x_len, 5*y_len))


        i = 0
        iterator = list(self.grids.items())
        for a in range(x_len):
            for b in range(y_len):
                if i == len(self.grids):
                    im = axes[a, b].imshow(self.temperature_grid, cmap='hot', interpolation='nearest')
                    axes[a, b].set_title(f'Temperature (t={self.time_step})')
                    axes[a, b].set_xlabel('X Position')
                    axes[a, b].set_ylabel('Y Position')
                    plt.colorbar(im, ax=axes[a, b], label=f'Temperature [K]')
                    break
                else:
                    source_type, grid = iterator[i]
                    i += 1
                    im = axes[a, b].imshow(grid, cmap='hot', interpolation='nearest', vmin=0)
                    axes[a, b].set_title(f'{source_type.title()} Distribution (t={self.time_step})')
                    axes[a, b].set_xlabel('X Position')
                    axes[a, b].set_ylabel('Y Position')
                    plt.colorbar(im, ax=axes[a, b], label=f'{source_type} amount [mol]', use_gridspec=True)

                    for source in [src for src in self.sources if src.type == source_type]:
                        axes[a, b].plot(source.x, source.y, 'x', markersize=8)

        return fig, axes


if __name__ == "__main__":
    env = Environment(width=16, height=16)
    env.add_source(x=5, y=5, intensity=10.0, radius=3)
    env.add_source(x=5, y=5, intensity=0.1, radius=3, resource_type='C1O2')
    env.add_source(x=10, y=10, intensity=0.1, radius=1, resource_type='C1O2')

    for x in range(500):
        env.step()
        compounds_extracted = env.extract_compounds_from(5, 5, 1.0)
        energy_extracted = env.extract_energy_from(5, 5, 0.7)
        if x % 100 == 0:
            fig = env.visualize()
            plt.show()
    fig = env.visualize()
    plt.show()