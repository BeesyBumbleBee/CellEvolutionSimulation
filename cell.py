from __future__ import annotations
from copy import deepcopy

import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass
from chemistry import Atom, Compound, ReactionEngine
from environment import Environment


class PossibleReactant:
    def __init__(self,
                 need_atom: str|None = None):
        self.need_atom = need_atom

    def check_reactant(self, reactant: str) -> bool:
        if self.need_atom:
            return self.need_atom in reactant
        return False

    def __repr__(self):
        return f"A compound with '{self.need_atom}' atom"


class ReactionRule:
    def __init__(self, num_reactants: int, possible_reactants: List[PossibleReactant], use_energy: int, priority: float = 1.0):
        assert num_reactants > 0
        assert len(possible_reactants) == num_reactants
        assert use_energy >= 0

        self.num_reactants = num_reactants
        self.possible_reactants = possible_reactants
        self.use_energy = use_energy

        self.priority = priority

    def __repr__(self):
        out_str = f"({self.priority:4.2f})Use {self.use_energy:6.2f} kJ/mol to connect ({self.num_reactants}) reactants:\n"
        for reactant in self.possible_reactants:
            out_str += f'\t{reactant}\n'
        return out_str

class Genome:
    def __init__(self, rules: List[ReactionRule]):
        self.rules = rules

    def mutate(self, mutation_rate: float = 0.4, rng: np.random = np.random.default_rng()) -> Genome:
        from copy import deepcopy
        possible_mutations = ['delete', 'change_num', 'change_energy', 'change_possible', 'change_priority']

        new_rules = []
        for rule in self.rules:
            new_rule = deepcopy(rule)
            if mutation_rate < rng.random():
                new_rules.append(new_rule)
                continue

            mutation = rng.choice(possible_mutations)
            if mutation == 'delete':
                continue

            elif mutation == 'change_num':
                new_num_reactants = rng.choice([x for x in range(2, 4) if x != new_rule.num_reactants])
                if new_num_reactants < new_rule.num_reactants:
                    new_rule.possible_reactants = new_rule.possible_reactants[:new_num_reactants]
                else:
                    new_rule.possible_reactants.extend([
                        PossibleReactant(rng.choice(list(Atom.basic_atoms.keys())))
                        for _ in range(new_num_reactants - new_rule.num_reactants)
                    ])
                new_rule.num_reactants = new_num_reactants

            elif mutation == 'change_energy':
                delta_energy = rng.uniform(10, 100)
                if rng.random() > 0.5:
                    new_rule.use_energy += delta_energy
                else:
                    new_rule.use_energy -= delta_energy
                if new_rule.use_energy < 0:
                    new_rule.use_energy = 0

            elif mutation == 'change_possible':
                change_idx = rng.choice([i for i in range(len(new_rule.possible_reactants))])
                new_rule.possible_reactants[change_idx] = PossibleReactant(rng.choice(list([x for x in Atom.basic_atoms.keys() if x != new_rule.possible_reactants[change_idx].need_atom])))

            elif mutation == 'change_priority':
                delta_priority = rng.uniform(0.1, 1.0)
                if rng.random() > 0.5:
                    new_rule.priority += delta_priority
                else:
                    new_rule.priority -= delta_priority
                if new_rule.priority < 0:
                    new_rule.priority = 0

            new_rules.append(new_rule)

        return Genome(new_rules)


class Protocell:
    next_id: int = 0

    def __init__(self, x: int, y: int, genome: Genome, initial_energy: float = 1000.0):
        self.id = Protocell.next_id
        Protocell.next_id += 1

        self.x = x
        self.y = y
        self.genome = genome
        self.energy = initial_energy

        self.compounds: Dict[str, List[Compound | float]] = {}

        self.age = 0
        self.alive = True
        self.generation = 0

        self.reproduction_cost = 800.0
        self.reproduction_threshold = np.random.uniform(self.reproduction_cost*1.5, self.reproduction_cost*3)
        self.base_metabolism = 5.0

    def __repr__(self):
        return f'Cell-{self.id:<4d} (x={self.x:<3d}, y={self.y:<3d}): Age: {self.age:4d} | Gen: {self.generation:2d} | Energy: {self.energy:8.2f} | ALIVE: {self.alive}'

    def summary(self) -> str:
        out_str = f"{self}\n"
        out_str += "Compounds:\n"
        for formula, val in self.compounds.items():
            out_str += f'\t{formula}: {val[1]:4.2f}\n'

        out_str += "Genome:\n"
        for rule in self.genome.rules:
            out_str += f"\t{rule}"
        return out_str



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
        self.add_compounds(absorbed)

    def add_compounds(self, compounds: dict[str, List[Compound | float]]):
        for compound, val in compounds.items():
            if compound not in self.compounds.keys():
                self.compounds[compound] = val
            else:
                self.compounds[compound][1] = val[1]

    def get_reactants(self, rule: ReactionRule):
        compounds = [
            [x for x in self.compounds if possible.check_reactant(x)]
            for possible in rule.possible_reactants
        ]
        return [np.random.choice(list(x)) if len(x) != 0 else None for x in compounds]


    def execute_reaction_rule(self, rule: ReactionRule):
        from copy import deepcopy

        reactants = self.get_reactants(rule)
        provided_energy = rule.use_energy if rule.use_energy < self.energy else self.energy
        used_reactants = {reactant_formula : round(self.compounds[reactant_formula][1]) for reactant_formula in reactants}

        reaction = ReactionEngine(provided_energy, energy_loss=0.2)

        for reactant_formula, count in used_reactants.items():
            reaction.add_reactants([deepcopy(self.compounds[reactant_formula][0]) for _ in range(count)])

        products, energy_released = reaction.evaluate_reaction()
        self.energy = energy_released - provided_energy
        count_products = {product.symbol : 0 for product in products}
        for product in products:
            count_products[product.symbol] += 1

        products = {product.symbol: [product, count_products[product.symbol]] for product in products}
        self.add_compounds(products)

    def execute_metabolism(self):
        for rule in self.genome.rules:
            self.execute_reaction_rule(rule)
        return

    def metabolise(self):
        self.energy -= self.base_metabolism
        if self.energy <= 0:
            self.energy = 0
            self.alive = False


