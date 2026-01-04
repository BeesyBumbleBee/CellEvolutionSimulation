from __future__ import annotations
from copy import deepcopy

import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass
from chemistry import Atom, Compound, ReactionEngine, CompoundStability
from environment import Environment


class PossibleReactant:
    def __init__(self, need_atoms: List[str] = None):
        self.need_atoms: List[str] = need_atoms
        if need_atoms is None:
            self.need_atoms = []


    def check_reactant(self, reactant: str) -> bool:
        if self.need_atoms:
            reactant = list(reactant)
            for atom in self.need_atoms:
                if atom not in reactant:
                    return False
                reactant.remove(atom)

        return True

    def __repr__(self):
        out_str = "Compound "
        if self.need_atoms:
            out_str += "with "
            for i, atom in enumerate(self.need_atoms):
                out_str += f"{atom} "
                if i < len(self.need_atoms) - 1:
                    out_str += "and "

        return out_str

    def get_mutated(self, mutation_rate: float = 1.0, rng: np.random = np.random.default_rng()) -> PossibleReactant:
        possible_mutations = ['add needed atom', 'remove needed atom']
        new_need_atoms = self.need_atoms

        if mutation_rate > rng.random():
            mutation = rng.choice(possible_mutations)

            if mutation == 'add needed atom':
                new_need_atoms.append(rng.choice(list(Atom.basic_atoms.keys())))
            elif mutation == 'remove needed atom':
                if len(self.need_atoms) != 0:
                    new_need_atoms.pop(rng.integers(len(self.need_atoms)))

        return PossibleReactant(
            need_atoms=new_need_atoms,
        )

    @staticmethod
    def get_random(rng: np.random = np.random.default_rng()) -> PossibleReactant:
        num_need_atoms = rng.integers(1, 3)
        return PossibleReactant([rng.choice(list(Atom.basic_atoms.keys())) for _ in range(num_need_atoms)])



class ReactionRule:
    def __init__(self, possible_reactants: List[PossibleReactant], use_energy: int, priority: float = 1.0):
        assert use_energy >= 0

        self.possible_reactants = possible_reactants
        self.use_energy = use_energy

        self.priority = priority

    def __repr__(self):
        out_str = f"Use {self.use_energy:6.2f} kJ/mol to connect ({len(self.possible_reactants)}) reactants:\n"
        for reactant in self.possible_reactants:
            out_str += f'\t{reactant}\n'
        return out_str

    def get_mutated(self, mutation_rate: float = 1.0, rng: np.random = np.random.default_rng()):
        possible_mutations = ['add reactant', 'remove reactant', 'mutate possible reactant', 'change used energy', 'change priority']

        new_possible_reactants = self.possible_reactants
        new_use_energy = self.use_energy
        new_priority = self.priority

        if mutation_rate > rng.random():
            mutation = rng.choice(possible_mutations)

            if mutation == 'add reactants':
                new_possible_reactants.append(PossibleReactant.get_random(rng=rng))

            elif mutation == 'remove reactants':
                if len(new_possible_reactants) > 0:
                    new_possible_reactants.pop(rng.integers(len(new_possible_reactants)))

            elif mutation == 'mutate possible reactant':
                if len(new_possible_reactants) > 0:
                    change_idx = rng.integers(len(new_possible_reactants))
                    new_possible_reactants[change_idx] = new_possible_reactants[change_idx].get_mutated(rng=rng)

            elif mutation == 'change used energy':
                if rng.random() > 0.5:
                    new_use_energy -= rng.integers(100)
                    if new_use_energy < 0:
                        new_use_energy = 0
                else:
                    new_use_energy += rng.integers(100)

            elif mutation == 'change priority':
                if rng.random() > 0.5:
                    new_priority -= rng.uniform(0.5)
                    if new_priority < 0:
                        new_priority = 0
                else:
                    new_priority += rng.uniform(0.5)


        return ReactionRule(
            possible_reactants=new_possible_reactants,
            use_energy=new_use_energy,
            priority=new_priority,
        )

    @staticmethod
    def get_random(rng: np.random = np.random.default_rng()) -> ReactionRule:
        num_reactants = rng.integers(2, 5)
        return ReactionRule(
            possible_reactants=[PossibleReactant.get_random(rng)
                                for _ in range(num_reactants)],
            use_energy=int(rng.integers(0, 300)),
            priority=round(rng.uniform(low=0.5, high=2.0), 2))


class Genome:
    def __init__(self, rules: List[ReactionRule]):
        self.rules = rules


    def get_mutated(self, mutation_rate: float = 0.4, rng: np.random = np.random.default_rng()) -> Genome:
        possible_mutations = ['add reaction rule', 'remove reaction rule', 'mutate reaction rule']
        new_rules = self.rules

        if mutation_rate > rng.random():
            mutation = rng.choice(possible_mutations)

            if mutation == 'add reaction rule':
                new_rules.append(ReactionRule.get_random(rng=rng))
            elif mutation == 'remove reaction rule':
                if len(new_rules) > 0:
                    new_rules.pop(rng.integers(len(new_rules)))
            elif mutation == 'mutate reaction rule':
                if len(new_rules) > 0:
                    change_idx = rng.integers(len(new_rules))
                    new_rules[change_idx] = new_rules[change_idx].get_mutated(rng=rng)
        return Genome(
            rules=new_rules
        )

    @staticmethod
    def get_random(rng: np.random = np.random.default_rng()) -> Genome:
        num_rules = rng.integers(1, 2)
        return Genome(rules=[ReactionRule.get_random(rng) for _ in range(num_rules)])


class Protocell:
    next_id: int = 0
    rng: np.random = np.random.default_rng()

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
        self.reproduction_threshold = Protocell.rng.normal(self.reproduction_cost*2, self.reproduction_cost*0.5)
        self.base_metabolism = 10.0

        self.total_decompositions = 0
        self.energy_from_environment = 0.0
        self.energy_from_reactions = 0.0
        self.energy_from_decomposition = 0.0

    def __repr__(self):
        return f'Cell-{self.id:<4d} (x={self.x:<3d}, y={self.y:<3d}): Age: {self.age:4d} | Gen: {self.generation:2d} | Energy: {self.energy:8.2f} | ALIVE: {self.alive}'

    def summary(self) -> str:
        out_str = f"{self}\n"
        out_str += f"Energy from environment: {self.energy_from_environment:8.2f}\n"
        out_str += f"Energy from reactions: {self.energy_from_reactions:8.2f}\n"
        out_str += f"Energy from decomposition: {self.energy_from_decomposition:8.2f}\n"
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

        energy = env.extract_energy_from(self.x, self.y, absorption_rate)
        self.energy_from_environment += energy
        self.energy += energy

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
        return [x for x in (Protocell.rng.choice(list(x)) if len(x) != 0 else None for x in compounds) if x]

    def check_compound_stability(self, instability_threshold: float = 0.5):
        """
        Check all compounds for instability and decompose if needed.

        This is called after metabolism to prevent runaway chain formation.
        """
        new_compounds = {}
        total_energy_released = 0.0
        decomposition_count = 0

        for formula, (compound, count) in list(self.compounds.items()):
            instability = CompoundStability.calculate_instability(compound)

            # Decompose if too unstable
            if instability > instability_threshold and not compound.stable:
                decomposition_count += int(count)

                # Decompose each instance
                for _ in range(int(count)):
                    fragments, energy = CompoundStability.decompose_compound(
                        deepcopy(compound),
                        max_breaks=min(3, int(instability) + 1)
                    )
                    total_energy_released += energy

                    # Add fragments
                    for fragment in fragments:
                        frag_formula = fragment.symbol
                        if frag_formula not in new_compounds:
                            new_compounds[frag_formula] = [fragment, 0]
                        new_compounds[frag_formula][1] += 1
            else:
                # Keep stable compound
                if formula not in new_compounds:
                    new_compounds[formula] = [compound, 0]
                new_compounds[formula][1] += count

        # Update compounds and energy
        self.compounds = new_compounds
        self.energy += total_energy_released

        self.total_decompositions += decomposition_count
        self.energy_from_decomposition += total_energy_released

        return decomposition_count, total_energy_released

    def execute_reaction_rule(self, rule: ReactionRule):
        from copy import deepcopy

        reactants = self.get_reactants(rule)
        provided_energy = rule.use_energy if rule.use_energy < self.energy else self.energy
        used_reactants = {reactant_formula : min(round(self.compounds[reactant_formula][1]), 2) for reactant_formula in reactants}

        reaction = ReactionEngine(provided_energy, energy_loss=0.2)

        for reactant_formula, count in used_reactants.items():
            reaction.add_reactants([deepcopy(self.compounds[reactant_formula][0]) for _ in range(count)])

        products, energy_released = reaction.evaluate_reaction()
        self.energy += energy_released - provided_energy
        self.energy_from_reactions += energy_released - provided_energy
        count_products = {product.symbol : 0 for product in products}
        for product in products:
            count_products[product.symbol] += 1

        products = {product.symbol: [product, count_products[product.symbol]] for product in products}
        self.add_compounds(products)

    def execute_metabolism(self):
        for rule in self.genome.rules:
            self.execute_reaction_rule(rule)

        decompositions, energy_recovered = self.check_compound_stability(
            instability_threshold=0.8
        )

        return

    def metabolise(self):
        self.energy -= self.base_metabolism * max((self.age / 32), 1)
        if self.energy <= 0:
            self.energy = 0
            self.alive = False


