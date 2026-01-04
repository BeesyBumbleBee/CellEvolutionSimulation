import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass
from chemistry import Atom, Compound, ReactionEngine
from environment import Environment

class ReactionRule:
    def __init__(self, num_reactants: int, use_energy: int):
        pass


class Genome:
    def __init__(self, rules: List):
        pass

class Protocell:
    next_id: int = 0

    def __init__(self, x: int, y: int, genome: Genome, initial_energy: float = 1000.0):
        self.id = Protocell.next_id
        Protocell.next_id += 1

        self.x = x
        self.y = y
        self.genome = genome
        self.energy = initial_energy

        self.compounds = {}

        self.age = 0
        self.alive = True
        self.generation = 0

        self.reproduction_cost = 800.0
        self.reproduction_threshold = np.random.uniform(self.reproduction_cost*1.5, self.reproduction_cost*3)
        self.base_metabolism = 5.0

    def __repr__(self):
        return f'Cell-{self.id:<4d} (x={self.x:<3d}, y={self.y:<3d}): Age: {self.age:4d} | Gen: {self.generation:2d} | Energy: {self.energy:8.2f} | ALIVE: {self.alive}'

    @property
    def can_reproduce(self):
        return self.energy >= self.reproduction_threshold


    def step(self, env: Environment):
        if not self.alive:
            return

        self.age += 1
        # 1. Absorb energy from environment
        self.absorb_energy(env)

        # 2. Absorb compounds from environment
        self.absorb_compounds(env)

        # 3. Execute reactions listed in genome
        self.execute_metabolism()

        # 4. Metabolise -> use energy to keep alive
        self.metabolise()

        return

    def absorb_energy(self, env: Environment, absorption_rate: float = 0.7):
        if absorption_rate > 1.0:
            absorption_rate = 1.0
        if absorption_rate < 0.0:
            absorption_rate = 0.0

        self.energy += env.extract_energy_from(self.x, self.y, absorption_rate)

    def absorb_compounds(self, env: Environment, absorption_rate: float = 0.7):
        if absorption_rate > 1.0:
            absorption_rate = 1.0
        if absorption_rate < 0.0:
            absorption_rate = 0.0

        absorbed = env.extract_compounds_from(self.x, self.y, absorption_rate)
        for compound, val in absorbed.items():
            if compound not in self.compounds.keys():
                self.compounds[compound] = val
            else:
                self.compounds[compound] += val

    def execute_metabolism(self):
        return

    def metabolise(self):
        self.energy -= self.base_metabolism
        if self.energy <= 0:
            self.energy = 0
            self.alive = False


