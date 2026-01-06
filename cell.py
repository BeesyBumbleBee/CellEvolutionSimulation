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


    def check_reactant(self, reactant: Compound) -> bool:
        if self.need_atoms:
            reactant = list(reactant.symbol)
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
        from copy import deepcopy
        possible_mutations = ['add needed atom', 'remove needed atom']
        new_need_atoms = deepcopy(self.need_atoms)

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
        return f"Reaction rule: Energy to use = {self.use_energy}, Possible reactants = {self.possible_reactants}"

    def summary(self) -> str:
        out_str = f"Use {self.use_energy:6.2f} kJ/mol to connect ({len(self.possible_reactants)}) reactants:\n"
        for reactant in self.possible_reactants:
            out_str += f'\t{reactant}\n'
        return out_str

    def get_mutated(self, mutation_rate: float = 1.0, rng: np.random = np.random.default_rng()):
        from copy import copy
        possible_mutations = ['add reactant', 'remove reactant', 'mutate possible reactant', 'change used energy', 'change priority']

        new_possible_reactants = copy(self.possible_reactants)
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
        num_reactants = rng.integers(2, 4)
        return ReactionRule(
            possible_reactants=[PossibleReactant.get_random(rng)
                                for _ in range(num_reactants)],
            use_energy=int(rng.integers(0, 300)),
            priority=round(rng.uniform(low=0.5, high=2.0), 2))


class Genome:
    def __init__(self, rules: List[ReactionRule]):
        self.rules = rules


    def get_mutated(self, mutation_rate: float = 0.4, rng: np.random = np.random.default_rng()) -> Genome:
        from copy import deepcopy
        possible_mutations = ['add reaction rule', 'remove reaction rule', 'mutate reaction rule']
        new_rules = deepcopy(self.rules)

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
        num_rules = 1
        return Genome(rules=[ReactionRule.get_random(rng) for _ in range(num_rules)])

@dataclass
class CellLog:
    age: int = -1
    generation: int = -1
    energy_used: int = 0
    energy_from_reaction: int = 0
    energy_from_decomposition: int = 0
    energy_from_environment: int = 0
    number_of_reactions: int = 0
    number_of_decompositions: int = 0
    absorbed_amount: int = 0
    reactants_used: List[Tuple[str]] = None
    reactions: List[str] = None


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

        self.compounds: List[Compound] = []

        self.age = 0
        self.alive = True
        self.generation = 0

        self.reproduction_cost = 800.0
        self.reproduction_threshold = Protocell.rng.normal(self.reproduction_cost*2, self.reproduction_cost*0.1)
        self.base_metabolism = 5.0

        self.total_decompositions = 0
        self.energy_from_environment = 0.0
        self.energy_from_reactions = 0.0
        self.energy_from_decomposition = 0.0
        self.reactions_count: Dict[str, int] = {}
        self.children = 0
        self.parent: str = ""

    @property
    def unique_compounds(self) -> Dict[str, Tuple[Compound, int]]:
        return {
            compound.symbol:
                (
                    compound,
                    len([x for x in self.compounds if x.__hash__() == compound.__hash__()])
                )
                for compound in self.compounds
        }

    def __repr__(self):
        return f'Cell-{self.id:<4d} (x={self.x:<3d}, y={self.y:<3d}): Age: {self.age:4d} | Gen: {self.generation:2d} | Energy: {self.energy:8.2f} | ALIVE: {self.alive}'

    def summary(self) -> str:
        out_str = f"{self}\n"
        out_str += f'Parent: {self.parent:<10s} | Children: {self.children:<4d}\n'
        out_str += f"Energy from environment:    {self.energy_from_environment:>8.2f}\n"
        out_str += f"Energy from reactions:      {self.energy_from_reactions:>8.2f}\n"
        out_str += f"Energy from decomposition:  {self.energy_from_decomposition:>8.2f}\n"
        out_str += f"Most common reactions: \n"
        sorted_reaction_counts = list(sorted([(key, val) for key, val in self.reactions_count.items()], key=lambda x: x[1], reverse=True))
        for i, x in enumerate(sorted_reaction_counts[:5 if len(sorted_reaction_counts) > 5 else len(sorted_reaction_counts)]):
            reaction, val = x
            out_str += f'\t[{i:>2d}] ({val: 4d}) {reaction}'

        out_str += "\nCompounds:\n"
        for compound_formula, val in self.unique_compounds.items():
            out_str += f'\t{compound_formula:>10s} : {val[1]:3d} | {val[0]}\n'

        out_str += "\nGenome:\n"
        for i, rule in enumerate(self.genome.rules):
            out_str += f"[{i:>2d}] {rule.summary()}\n"
        return out_str



    @property
    def can_reproduce(self):
        return self.energy >= self.reproduction_threshold


    def step(self, env: Environment) -> CellLog:
        log = {}
        if not self.alive:
            return CellLog()

        self.age += 1
        log['age'] = self.age
        log['generation'] = self.generation

        # 1. Absorb energy from environment
        log['energy_from_environment'] = self.absorb_energy(env)

        # 2. Absorb compounds from environment
        log['absorbed_amount'] = self.absorb_compounds(env)

        # 3. Execute reactions listed in genome
        log.update(self.execute_metabolism())

        # 4. Metabolise -> use energy to keep alive
        log['energy_used'] += self.metabolise()

        return CellLog(
            **log
        )

    def absorb_energy(self, env: Environment, absorption_rate: float = 0.7) -> float:
        if absorption_rate > 1.0:
            absorption_rate = 1.0
        if absorption_rate < 0.0:
            absorption_rate = 0.0

        energy = env.extract_energy_from(self.x, self.y, absorption_rate)
        self.energy_from_environment += energy
        self.energy += energy
        return energy

    def absorb_compounds(self, env: Environment, absorption_rate: float = 1.0):
        if absorption_rate > 1.0:
            absorption_rate = 1.0
        if absorption_rate < 0.0:
            absorption_rate = 0.0

        absorbed = env.extract_compounds_from(
            self.x,
            self.y,
            absorption_rate,
            {x.symbol: len([y for y in self.compounds if y.symbol == x.symbol]) for x in self.compounds} # only absorb if have less of compound than in environment
        )
        self.compounds.extend(absorbed)
        return len(absorbed)

    def get_reactants(self, rule: ReactionRule):
        compounds = []
        for possible_reactant in rule.possible_reactants:
            possible_compounds = list([x for x in self.compounds if possible_reactant.check_reactant(x)])
            added_compound = Protocell.rng.choice(possible_compounds) if possible_compounds else None
            if added_compound:
                self.compounds.remove(added_compound)
                compounds.append(added_compound)
        return compounds

    def check_compound_stability(self, instability_threshold: float = 0.5):
        new_compounds = []
        total_energy_released = 0.0
        decomposition_count = 0

        new_compounds.extend(list([x for x in self.compounds if x.stable]))

        for compound in [x for x in self.compounds if not x.stable]:
            instability = CompoundStability.calculate_instability(compound)

            # Decompose if too unstable
            if instability > instability_threshold:
                decomposition_count += 1

                fragments, energy = CompoundStability.decompose_compound(
                    deepcopy(compound),
                    max_breaks=min(3, int(instability) + 1)
                )
                total_energy_released += energy

                new_compounds.extend(fragments)
            else:
                new_compounds.append(compound)

        self.compounds = new_compounds
        self.energy += total_energy_released

        self.total_decompositions += decomposition_count
        self.energy_from_decomposition += total_energy_released

        return decomposition_count, total_energy_released

    def execute_reaction_rule(self, rule: ReactionRule) -> Tuple[int, int, Tuple[str], str] | Tuple[int, int, None, None]:
        reactants = self.get_reactants(rule)
        if len(reactants) <= 1:
            return 0, 0, None, None
        provided_energy = rule.use_energy if rule.use_energy < self.energy else self.energy

        reaction = ReactionEngine(provided_energy, energy_loss=0.2)
        reaction.add_reactants(reactants)

        products, energy_released = reaction.evaluate_reaction()
        self.energy += energy_released - provided_energy
        self.energy_from_reactions += energy_released - provided_energy

        self.compounds.extend(products)

        return provided_energy, energy_released, tuple((x.symbol for x in reactants)), reaction.reaction_summary

    def execute_metabolism(self):
        total_energy_used = 0
        total_energy_from_reaction = 0
        all_reactants = []
        reactions = []
        number_of_reactions = 0
        for rule in self.genome.rules:
            energy_used, energy_from_reaction, reactants_used, reaction_summary = self.execute_reaction_rule(rule)
            if reactants_used is not None:
                total_energy_used += energy_used
                total_energy_from_reaction += energy_from_reaction
                all_reactants.append(reactants_used)
                reactions.append(reaction_summary)
                try:
                    self.reactions_count[reaction_summary] += 1
                except KeyError:
                    self.reactions_count[reaction_summary] = 1
                number_of_reactions += 1

        decompositions, energy_recovered = self.check_compound_stability(
            instability_threshold=0.8
        )

        return {
            'number_of_decompositions': decompositions,
            'energy_from_decomposition': energy_recovered,
            'energy_used': total_energy_used,
            'number_of_reactions': number_of_reactions,
            'energy_from_reaction': total_energy_from_reaction,
            'reactants_used': all_reactants,
            'reactions': reactions,
        }

    def metabolise(self):
        self.energy -= self.base_metabolism * max((self.age / 32), 1) # apply pressure to reproduce and not stagnate
        if self.energy <= 0:
            self.energy = 0
            self.alive = False

        return self.base_metabolism * max((self.age / 32), 1)


