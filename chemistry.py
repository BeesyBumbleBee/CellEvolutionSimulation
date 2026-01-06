from __future__ import annotations
from enum import StrEnum
from itertools import combinations
from typing import List, Dict, Optional, Tuple, Generator
from dataclasses import dataclass
import logging
import networkx as nx

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


class CompoundStability:
    """
    Tracks and manages compound instability.

    Instability factors:
    1. Unsatisfied electrons (radicals)
    2. Extreme chain length
    3. Strained bonds
    4. High partial charges
    """

    @staticmethod
    def calculate_instability(compound) -> float:
        """Calculate overall instability score (0.0 = stable, 1.0+ = very unstable)"""

        # Factor 1: Unsatisfied electrons (most important)
        radical_penalty = sum(1 for atom in compound.components if atom.electrons_needed > 0) * 0.3

        # Factor 2: Chain length penalty (exponential for very long chains)
        chain_length = len(compound.components)
        if chain_length > 10:
            length_penalty = (chain_length - 10) * 0.05
        else:
            length_penalty = 0.0

        # Factor 3: Bond strain
        strain_penalty = sum(bond.strain for bond, _, _ in compound.bonds) * 0.2

        # Factor 4: Extreme partial charges
        charge_penalty = sum(abs(atom.partial_charge) for atom in compound.components
                             if abs(atom.partial_charge) > 1.0) * 0.1

        total_instability = radical_penalty + length_penalty + strain_penalty + charge_penalty

        return total_instability

    @staticmethod
    def should_decompose(compound, instability_threshold: float = 0.5) -> bool:
        """Determine if compound should spontaneously decompose"""
        instability = CompoundStability.calculate_instability(compound)
        return instability > instability_threshold

    @staticmethod
    def find_weakest_bonds(compound, n: int = 3) -> List[Tuple[int, float]]:
        """Find indices of n weakest bonds for decomposition"""
        if len(compound.bonds) == 0:
            return []

        bond_strengths = [
            (idx, bond.effective_energy)
            for idx, (bond, _, _) in enumerate(compound.bonds)
        ]

        sorted_bonds = sorted(bond_strengths, key=lambda x: x[1])
        return sorted_bonds[:min(n, len(sorted_bonds))]

    @staticmethod
    def make_inner_connections(compound: Compound) -> Tuple[int, Compound]:
        reaction = ReactionEngine(0)
        reaction.add_reactant(compound)
        reaction.make_inner_reactant_connections()
        energy_gained = reaction.dissipate_energy()
        return energy_gained, reaction.reactants[0]

    @staticmethod
    def decompose_compound(compound, max_breaks: int = 2) -> Tuple[List, float]:
        if len(compound.bonds) == 0:
            return [compound], 0.0

        energy_gained, compound = CompoundStability.make_inner_connections(compound)
        if not CompoundStability.should_decompose(compound, instability_threshold=0.8):
            return [compound], energy_gained

        instability = CompoundStability.calculate_instability(compound)

        # Determine how many bonds to break based on instability
        num_breaks = min(max_breaks, int(instability) + 1)
        weak_bonds = CompoundStability.find_weakest_bonds(compound, num_breaks)

        if not weak_bonds:
            return [compound], 0.0

        energy_released = 0.0

        # Break bonds from highest index to lowest to avoid index shifting
        for bond_idx, bond_energy in sorted(weak_bonds, reverse=True):
            bond, i, j = compound.bonds[bond_idx]

            # Release partial energy (unstable decomposition is less efficient)
            energy_released += bond_energy * 0.1

            # Remove bond
            compound.bonds.pop(bond_idx)

            # Update electron counts
            bond.component_A.cov_electrons -= bond.multiplicity
            bond.component_B.cov_electrons -= bond.multiplicity

        # Split into separate compounds
        fragments = compound.split_compounds()

        # Distribute remaining energy among fragments
        if compound.remaining_energy > 0:
            total_mass = sum(f.mass for f in fragments)
            for fragment in fragments:
                fragment.remaining_energy = int(
                    (fragment.mass / total_mass) * compound.remaining_energy
                )

        return fragments, energy_released


class Atom:
    class AtomNotDefined(Exception):
        pass

    atoms_electronegativity = {
        'H': 2.20,
        'C': 2.55,
        'N': 3.04,
        'O': 3.44,
    }

    basic_atoms: Dict[str, Dict] = {
        'H': {'mass': 2, 'symbol': 'H', 'electrons_in_covalence': 1, 'optimal_electrons': 2},
        'C': {'mass': 12, 'symbol': 'C', 'electrons_in_covalence': 4, 'optimal_electrons': 8},
        'N': {'mass': 14, 'symbol': 'N', 'electrons_in_covalence': 5, 'optimal_electrons': 8},
        'O': {'mass': 16, 'symbol': 'O', 'electrons_in_covalence': 6, 'optimal_electrons': 8},
    }

    def __init__(self, mass: float, symbol: str, electrons_in_covalence: int, optimal_electrons: int):
        assert optimal_electrons >= electrons_in_covalence
        self.mass: float = mass
        self.symbol: str = symbol
        self.cov_electrons: int = electrons_in_covalence
        self.optimal_electrons: int = optimal_electrons
        self.partial_charge: float = 0.0

    def __repr__(self):
        status = "✓" if self.cov_electrons >= self.optimal_electrons else "✗"
        charge_str = f" δ= {self.partial_charge:+.2f}" if abs(self.partial_charge) > 0.01 else ""
        radical_str = "•" if self.is_radical else " "
        return f'{self.symbol:>2s}{radical_str} | Mass: {self.mass:5.2f} | Electrons: {self.cov_electrons:2d}/{self.optimal_electrons:2d} ({status}) |  Partial Charge: {charge_str}'

    @property
    def is_radical(self):
        return self.cov_electrons % 2 == 1

    @property
    def electronegativity(self) -> float:
        try:
            return Atom.atoms_electronegativity[self.symbol]
        except KeyError:
            raise Atom.AtomNotDefined

    @property
    def electrons_needed(self) -> int:
        return self.optimal_electrons - self.cov_electrons

    @staticmethod
    def get(atom_symbol: str) -> Atom:
        try:
            return Atom(**Atom.basic_atoms[atom_symbol.upper()])
        except KeyError:
            raise Atom.AtomNotDefined


class Bond:
    """
    Class defining energy values of bonds between two atoms
    Source: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Bond_Energies
    """
    class BondNotDefinedError(Exception):
        pass

    class AtomNotInBond(Exception):
        pass

    BondEnergy = {
        1: {
            'H-H': 432, 'H-O': 467, 'H-N': 391,
            'C-H': 413, 'C-C': 347, 'C-O': 358, 'C-N': 305,
            'O-O': 146,
            'N-N': 160, 'N-O': 201,

        },
        2: {
            'C-C': 614, 'C-O': 745, 'C-N': 615,
            'O-O': 495,
            'N-N': 418, 'N-O': 607,
        },
        3: {
            'C-C': 839, 'C-O': 1072, 'C-N': 891,
            'N-N': 941,
        }
    }

    def __init__(self, component_a: Atom, component_b: Atom, multiplicity: int = 1):
        self.component_A: Atom = component_a
        self.component_B: Atom = component_b
        if self.component_A.symbol > self.component_B.symbol:
            temp = self.component_A
            self.component_A = self.component_B
            self.component_B = temp
        self.multiplicity: int = multiplicity
        try:
            self.energy = Bond.BondEnergy[multiplicity][f'{self.component_A.symbol}-{self.component_B.symbol}']
        except KeyError:
            raise Bond.BondNotDefinedError
        self.strain: float = 0.0

    def __repr__(self):
        left_polar = '<' if self.polar_to_A else ' '
        right_polar = '>' if self.polar_to_B else ' '
        return f"{self.component_A.symbol:>2s}{left_polar}{['-', '=', '≡'][self.multiplicity - 1]}{right_polar}{self.component_B.symbol:2s} | Polarity {self.polarity:3.2f} | Energy: {self.energy} ({self.effective_energy})kJ/mol"

    @property
    def effective_energy(self) -> int:
        # TODO: Implement better function to get strained energy
        return round(self._strain_effect * self.energy)

    @property
    def _strain_effect(self) -> float:
        return 1.0 - (min(self.strain, 0.4) * min(self.polarity / 1.7, 1.2))

    @property
    def polarity(self) -> float:
        """
        Positive values mean polarity in direction of Atom A - negative Atom B
        """
        return abs(self.component_A.electronegativity - self.component_B.electronegativity)

    @property
    def polar_to_A(self) -> bool:
        return (self.component_A.electronegativity - self.component_B.electronegativity) > 0

    @property
    def polar_to_B(self) -> bool:
        return (self.component_B.electronegativity - self.component_A.electronegativity) > 0

    def comp_polarity(self, atom_symbol: str):
        if atom_symbol.upper() not in [self.component_A.symbol, self.component_B.symbol]:
            raise Bond.AtomNotInBond

        if atom_symbol == self.component_A.symbol:
            return -self.polarity if self.polar_to_A else self.polarity
        else:
            return -self.polarity if self.polar_to_B else self.polarity

    def decrease_mult(self) -> Tuple[int, bool]:
        """
        Decrease bond multiplity
        :return:
        energy_difference: int
        was_broken: bool
        """
        if self.multiplicity == 1:
            return self.effective_energy, True

        self.multiplicity -= 1
        old_energy = self.effective_energy
        self.energy = Bond.BondEnergy[self.multiplicity][f'{self.component_A.symbol}-{self.component_B.symbol}']
        return abs(old_energy-self.effective_energy), False

    def increase_mult(self) -> Tuple[int, bool]:
        """
        Increase bond multiplicity

        :return:
        energy_difference: int
        was_increased: bool
        """
        try:
            self.multiplicity += 1
            old_energy = self.effective_energy
            self.energy = Bond.BondEnergy[self.multiplicity][f'{self.component_A.symbol}-{self.component_B.symbol}']
            return abs(old_energy-self.effective_energy), True
        except KeyError:
            return 0, False

    @property
    def next_multiplicity_energy_difference(self) -> int:
        try:
            old_energy = self.effective_energy
            new_energy = Bond.BondEnergy[self.multiplicity+1][f'{self.component_A.symbol}-{self.component_B.symbol}']
            return abs(old_energy - round(new_energy*self._strain_effect))
        except KeyError:
            return -1

    @property
    def previous_multiplicity_energy_difference(self) -> int:
        if self.multiplicity == 1:
            return self.effective_energy
        old_energy = self.effective_energy
        new_energy = Bond.BondEnergy[self.multiplicity-1][f'{self.component_A.symbol}-{self.component_B.symbol}']
        return abs(old_energy - round(new_energy*self._strain_effect))

    @staticmethod
    def get_bond_energy(symbol_a: str, symbol_b: str, multiplicity: int) -> Optional[int]:
        if symbol_a > symbol_b:
            symbol_a, symbol_b = symbol_b, symbol_a
        try:
            return Bond.BondEnergy[multiplicity][f'{symbol_a}-{symbol_b}']
        except KeyError:
            return None


class Compound:
    def __init__(self, components: List[Atom],
                 preserve_bonds: Optional[List[Tuple[int, int, int]]] = None,
                 provided_energy: int = 0):
        """
        components: list of Component objects
        preserve_bonds: optional list of (idx_a, idx_b, multiplicity) bonds to preserve from reactants
        provided_energy: energy available for synthesis
        """
        self.components: List[Atom] = components
        self.bonds = []
        self.preserve_bonds = preserve_bonds or []
        self.remaining_energy: int = provided_energy

        self.mass = sum([x.mass for x in self.components])

        symbols = []
        symbols_count = {com[0]: len(com) for com in (list(x.symbol for x in self.components if x.symbol == y.symbol) for y in set(self.components))}
        list(map(lambda x: symbols.extend([str(x[0]), str(x[1])]), sorted(symbols_count.items(), key=lambda x: x[0], reverse=False)))
        self.symbol = "".join(symbols)

        self.optimal_electrons = sum([x.optimal_electrons for x in components])
        self.electrons = sum([x.cov_electrons for x in self.components])

    def __repr__(self):
        return f'{self.symbol:12s} | STABLE: {self.stable} | Mass: {self.mass:6.2f} | Energy remaining: {self.remaining_energy:8.2f} kJ/mol'

    @property
    def stable(self) -> bool:
        return all([x.electrons_needed == 0 for x in self.components])

    @property
    def graph(self) -> nx.Graph:
        graph = nx.Graph()
        graph.add_nodes_from([i for i in range(len(self.components))])
        graph.add_edges_from([[i, j] for _, i, j in self.bonds])
        return graph

    @property
    def connection_graphs(self) -> List[nx.Graph]:
        return [sub_graph for sub_graph in nx.connected_components(self.graph)]

    @property
    def is_connected(self) -> bool:
        return len(self.connection_graphs) == 1

    @property
    def reactive_sites(self) -> Dict[ReactiveSite.Character, List[ReactiveSite]]:
        return {
            ReactiveSite.Character.electrophilic: sorted([
                ReactiveSite(i, abs(atom.partial_charge), ReactiveSite.Character.electrophilic) for i, atom in
                enumerate(self.components) if atom.partial_charge > 0],
                key=lambda x: x.value
            ),
            ReactiveSite.Character.neutrophilic: sorted([
                ReactiveSite(i, abs(atom.partial_charge), ReactiveSite.Character.neutrophilic) for i, atom in
                enumerate(self.components) if atom.partial_charge < 0],
                key=lambda x: x.value
            ),
        }

    @property
    def most_reactive_site(self) -> ReactiveSite | None:
        if not self.is_reactive:
            return None
        max_val = 0
        most_reactive_site = None
        for character in ReactiveSite.Character:
            reactive_sites = self.reactive_sites[character]
            site = reactive_sites[0]
            if site.value > max_val:
                max_val = site.value
                most_reactive_site = site
        return most_reactive_site

    @property
    def is_reactive(self) -> bool:
        return all(self.reactive_sites.values())

    def __hash__(self):
        return hash(f'{self.symbol}' + '|'.join([f'{bond.component_A.symbol}{bond.multiplicity}{bond.component_B.symbol},{i},{j}' for bond,i,j in self.bonds]))

    def summary(self) -> str:
        out_str = '*' + '=' * 80 + '*\n'
        out_str += f"{self} \nStructure:\n"
        out_str += f"Atoms ({len(self.components)}):\n"
        for i, comp in enumerate(self.components):
            out_str += f"  [{i:2d}] {comp}\n"

        out_str += f"\nBonds ({len(self.bonds)}):\n"
        for i, bond in enumerate(self.bonds):
                out_str += f"  [{i:2d}] {bond[0]} kJ/mol)\n"
        out_str += '*' + '=' * 80 + '*\n'
        return out_str

    def draw_compound(self):
        import matplotlib.pyplot as plt


        single = list([(f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$', {'weight': bond.energy / 100}) for bond, i, j in self.bonds if bond.multiplicity >= 1])
        double = list([(f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$', {'weight': bond.energy / 100}) for bond, i, j in self.bonds if bond.multiplicity >= 2])
        triple = list([(f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$', {'weight': bond.energy / 100}) for bond, i, j in self.bonds if bond.multiplicity >= 3])

        nucleophiles = [f'${self.components[i].symbol}_'+'{'f'{i}'+'}$' for i in range(len(self.components)) if self.components[i].partial_charge < 0]
        electrophiles = [f'${self.components[i].symbol}_'+'{'f'{i}'+'}$' for i in range(len(self.components)) if self.components[i].partial_charge > 0]
        neutral = [f'${self.components[i].symbol}_'+'{'f'{i}'+'}$' for i in range(len(self.components)) if self.components[i].partial_charge == 0]

        graph = nx.Graph()
        graph.add_nodes_from([f'${self.components[i].symbol}_'+'{'f'{i}'+'}$' for i in range(len(self.components))])
        graph.add_edges_from(single)
        graph.add_edges_from(double)
        graph.add_edges_from(triple)
        pos = nx.spring_layout(graph, seed=42, k=2.0, iterations=1000)
        nx.draw_networkx_nodes(graph, pos, node_size=450, node_color='white', edgecolors="black", nodelist=neutral)
        nx.draw_networkx_nodes(graph, pos, node_size=450, node_color='white', edgecolors="red", nodelist=nucleophiles)
        nx.draw_networkx_nodes(graph, pos, node_size=450, node_color='white', edgecolors="blue", nodelist=electrophiles)
        nx.draw_networkx_labels(graph, pos, font_size=10)
        nx.draw_networkx_edges(graph, pos, arrows=True,edgelist=single, connectionstyle='arc3, rad = 0.1')
        nx.draw_networkx_edge_labels(graph, pos, edge_labels={tuple(edge[:2]): f'{i}' for i, edge in enumerate(single)},
                                     font_size=8, font_color='gray', rotate=False, bbox=dict(facecolor='white', alpha=0.3, edgecolor='white'),label_pos=0.5)

        nx.draw_networkx_edges(graph, pos,arrows=True, edgelist=double, connectionstyle='arc3, rad = 0.25')
        nx.draw_networkx_edges(graph, pos,arrows=True, edgelist=triple, connectionstyle='arc3, rad = 0.4')


        plt.show()

    def dissipate_energy(self, energy_percent: float = 1.0) -> int:
        energy_dissipated = round(self.remaining_energy * energy_percent)
        self.remaining_energy -= energy_dissipated
        return energy_dissipated

    def split_compounds(self) -> List[Compound]:
        compounds = []
        for sub_graph in self.connection_graphs:
            components = [self.components[i] for i in sub_graph]
            scale = {comp: i for comp, i in zip(sub_graph, range(len(sub_graph)))}
            bonds = [(bond,scale[i],scale[j]) for bond, i, j in self.bonds if i in sub_graph or j in sub_graph]
            energy = round(self.remaining_energy * (sum([self.components[comp].mass for comp in sub_graph]) / self.mass))  # divide energy based on mass

            compounds.append(Compound(components, provided_energy=energy))
            compounds[-1].bonds = bonds

        return compounds

    def update_partial_charges(self):
        for i, comp in enumerate(self.components):
            comp.partial_charge = sum([
                bond.comp_polarity(comp.symbol)
                for bond, a, b in self.bonds if a == i or b == i
            ])

    @staticmethod
    def atom(atom_symbol:str, provided_energy:int = 0) -> Compound:
        return Compound([Atom.get(atom_symbol)], provided_energy=provided_energy)

    @staticmethod
    def from_formula(formula: str, provided_energy: int=0) -> Compound:
        """
        Create a compound from a chemical formula string
        Example: from_formula("H2O", 2000) or from_formula("C6H12O6", 50000)
        """
        import re

        # Parse formula: C6H12O6 -> [('C', 6), ('H', 12), ('O', 6)]
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)

        compound = None
        reaction = ReactionEngine(provided_energy)
        for symbol, count in matches:
            count = int(count) if count else 1
            for _ in range(count):
                reaction.add_reactant(Compound.atom(symbol))
        product, _ = reaction.evaluate_reaction()

        return product[0]


@dataclass
class ReactiveSite:
    class Character(StrEnum):
        electrophilic = 'electrophilic'
        neutrophilic = 'neutrophilic'

    atom_index: int
    value: float
    character: ReactiveSite.Character

    @property
    def opposite_character(self) -> ReactiveSite.Character:
        if self.character == ReactiveSite.Character.electrophilic:
            return ReactiveSite.Character.neutrophilic
        else:
            return ReactiveSite.Character.electrophilic


@dataclass
class ReactionBridge:
    bond: Bond
    reactant_a_idx: int
    reactant_b_idx: int
    comp_a_idx: int
    comp_b_idx: int


class ReactionEngine:
    def __init__(self, starting_energy: int = 0, energy_loss: float = 0.2):
        self.system_energy: int = starting_energy
        self.energy_from_reactions: int = 0
        self.energy_loss: float = energy_loss
        self.reactants: List[Compound] = []
        self.reaction_summary: str = ""

    def __repr__(self):
        return f"Reaction with {len(self.reactants):2d} reactants. System energy {self.system_energy:8.2f}. Energy loss {self.energy_loss:4.2f}"

    @property
    def reactive_reactants(self) -> List[Compound]:
        return [reactant for reactant in self.reactants if reactant.is_reactive]

    @property
    def reactive_sites(self) -> List[Dict[ReactiveSite.Character, List[ReactiveSite]]]:
        return [comp.reactive_sites for comp in self.reactants]

    @property
    def summary(self) -> str:
        out_str: str = ""
        out_str += f"System energy: {self.system_energy:8.2f} kJ/mol\n"
        out_str += f"Reactants ({len(self.reactants)})\n"
        for i, reactant in enumerate(self.reactants):
            out_str += f"\t[{i:2d}] {reactant} ({len(reactant.bonds)} bonds):\n"
            for j, bond in enumerate(reactant.bonds):
                out_str += f"\t\t[{j:2d}] {bond[0]}\n"
        return out_str

    def __update_reaction_summary(self, reactants: List[str], products: List[str]) -> None:
        self.reaction_summary = ' + '.join(reactants) + ' -> ' + ' + '.join(products)

    def print_summary(self) -> None:
        print(self.summary)

    def add_reactant(self, reactant: Compound) -> None:
        self.reactants.append(reactant)
        self.system_energy += reactant.remaining_energy

    def add_reactants(self, reactant_list: List[Compound]) -> None:
        self.reactants.extend(reactant_list)
        self.system_energy += sum(x.remaining_energy for x in reactant_list)

    def evaluate_strains(self):
        for compound_a in self.reactive_reactants:
            site_a = compound_a.most_reactive_site

            for compound_b in self.reactive_reactants:
                if compound_a == compound_b:
                    continue
                site_b = compound_b.reactive_sites[site_a.opposite_character][0]
                strain = max(abs(site_a.value - site_b.value), 0)
                for bond in [bond for bond, i, j in compound_b.bonds if site_b.atom_index in [i,j]]:
                    bond.strain += strain

    def bond_breaker(self) -> Generator[int]:
        """
        Break specified bonds and return energy released
        bond_indices: list of indices in self.bonds to break
        Returns: energy released (positive value)
        """
        energy_consumed = 0
        bonds = []
        for reactant in self.reactants:
            bonds.extend([(bond_info, reactant, idx) for idx, bond_info in enumerate(reactant.bonds)])

        sorted_bond_energies = sorted(bonds, reverse=False,
                                      key=lambda x: x[0][0].previous_multiplicity_energy_difference)
        for comp_bond, compound, idx in sorted_bond_energies:
            bond, i, j = comp_bond

            yield bond.previous_multiplicity_energy_difference
            yield 0

            _, broke = bond.decrease_mult()
            compound.components[i].cov_electrons -= bond.multiplicity
            compound.components[j].cov_electrons -= bond.multiplicity
            if broke:
                compound.bonds.remove((bond, i, j))

    def break_bonds(self):
        bond_breaker = self.bond_breaker()
        for energy_required in bond_breaker:
            if (self.system_energy * 0.8) - energy_required <= 0:
                break
            bond_breaker.__next__()
            self.energy_from_reactions -= energy_required
            self.system_energy -= energy_required

        new_reactants = []
        [new_reactants.extend(y for y in x.split_compounds()) for x in self.reactants]
        self.reactants = new_reactants


    def find_best_bridge(self) -> ReactionBridge | None:
        unstable = [(i, x) for i, x in enumerate(self.reactants) if not x.stable]
        bridges = []
        for reactant_a, reactant_b in combinations(unstable, 2):
            idx_a, reactant_a = reactant_a
            idx_b, reactant_b = reactant_b

            unstable_a = [(i, idx_a, comp) for i, comp in enumerate(reactant_a.components) if
                          comp.electrons_needed != 0]
            unstable_b = [(i, idx_b, comp) for i, comp in enumerate(reactant_b.components) if
                          comp.electrons_needed != 0]

            for comp_a in unstable_a:
                for comp_b in unstable_b:
                    bond = Bond(comp_a[2], comp_b[2], 1)

                    bridges.append(ReactionBridge(bond, comp_a[1], comp_b[1], comp_a[0], comp_b[0]))
        bridges = sorted(bridges, reverse=True, key=lambda x: x.bond.energy)
        if len(bridges) == 0:
            return None
        else:
            return bridges[0]


    def connect_reactants(self, bridge: ReactionBridge) -> int:
        reactant_a = self.reactants[bridge.reactant_a_idx]
        reactant_b = self.reactants[bridge.reactant_b_idx]

        reactant_a.components[bridge.comp_a_idx].cov_electrons += 1
        reactant_b.components[bridge.comp_b_idx].cov_electrons += 1

        a_idx = bridge.comp_a_idx
        b_idx = len(reactant_a.components) + bridge.comp_b_idx

        comps = reactant_a.components
        comps.extend(reactant_b.components)

        bonds = [(bridge.bond, a_idx, b_idx)]
        bonds.extend([(bond, comps.index(bond.component_A) , comps.index(bond.component_B)) for bond,_,_ in reactant_a.bonds])
        bonds.extend([(bond, comps.index(bond.component_A) , comps.index(bond.component_B)) for bond,_,_ in reactant_b.bonds])

        energy = reactant_a.remaining_energy
        energy += reactant_b.remaining_energy
        energy += bridge.bond.energy

        new_reactant = Compound(comps)
        new_reactant.bonds = bonds
        new_reactant.energy = energy
        self.reactants.remove(reactant_a)
        self.reactants.remove(reactant_b)
        self.reactants.append(new_reactant)

        return bridge.bond.energy

    def dissipate_energy(self) -> int:
        total_mass = sum([x.mass for x in self.reactants])
        total_energy_dissipated = 0
        if self.energy_from_reactions > 0:
            for reactant in self.reactants:
                energy_gained = (reactant.mass / total_mass) * self.energy_from_reactions
                energy_dissipated = energy_gained * self.energy_loss
                total_energy_dissipated += round(energy_dissipated)
                reactant.remaining_energy += round(energy_dissipated)

        if self.energy_from_reactions < 0:
            overdue_energy = 0
            for reactant in self.reactants:
                energy_lost = round((reactant.mass / total_mass) * (-self.energy_from_reactions))
                if reactant.remaining_energy >= energy_lost + overdue_energy:
                    reactant.remaining_energy -= energy_lost + overdue_energy
                    overdue_energy = 0
                else:
                    overdue_energy += (energy_lost + overdue_energy) - reactant.remaining_energy
                    reactant.remaining_energy = 0

        return total_energy_dissipated


    def update_partial_charges(self):
        for reactant in self.reactants:
            reactant.update_partial_charges()

    def make_inner_reactant_connections(self):
        for reactant in self.reactants:
            unsatisfied_atoms = [(i,x) for i,x in enumerate(reactant.components) if x.electrons_needed != 0]
            while len(unsatisfied_atoms) > 1:
                possible_bonds = sorted([bond for bond, i, j in reactant.bonds
                                         if i in [x[0] for x in unsatisfied_atoms] and
                                         j in [x[0] for x in unsatisfied_atoms] and
                                         bond.next_multiplicity_energy_difference != -1],
                                        key=lambda x: x.next_multiplicity_energy_difference,
                                        reverse=True)
                if len(possible_bonds) == 0:
                    break
                energy_gained, ok = possible_bonds[0].increase_mult()
                if not ok:
                    break
                possible_bonds[0].component_A.cov_electrons += 1
                possible_bonds[0].component_B.cov_electrons += 1
                self.system_energy += energy_gained
                self.energy_from_reactions += energy_gained
                unsatisfied_atoms = [(i,x) for i,x in enumerate(reactant.components) if x.electrons_needed != 0]

    def remove_strain(self):
        for reactant in self.reactants:
            for bond, _, _ in reactant.bonds:
                bond.strain = 0

    def evaluate_reaction(self, verbose: bool=False) -> Tuple[List[Compound], int]:
        """
        1. Evaluate strain
        2. Break weak bonds until energy is used up
        3. Find possible bridges between reactants
        4. Connect reactants with best bridges
        5. Make inner reactant connections
        6. Dissipate energy -> loss
        7. Return list of products
        """
        logger.info("*===== Evaluating reaction =====*")

        logger.info("*== Reaction state at start ==*")
        logger.info(self.summary)
        starting_reactants = [x.symbol for x in self.reactants]

        logger.debug("-== Phase 1. Evaluating bond strains ==-")
        self.evaluate_strains()
        logger.debug(self.summary)

        logger.debug("-== Phase 2. Breaking bonds ==-")
        self.break_bonds()
        logger.debug(self.summary)

        logger.debug("-== Phase 3. Connecting reactants ==-")
        while (bridge := self.find_best_bridge()) is not None:
            energy_gained = self.connect_reactants(bridge)
            self.system_energy += energy_gained
            self.energy_from_reactions += energy_gained
            self.update_partial_charges()
        logger.debug(self.summary)

        logger.debug("-== Phase 4. Making bonds inside reactants (double, triple bonds) ==-")
        self.make_inner_reactant_connections()
        logger.debug(self.summary)

        if all(reactant.symbol in starting_reactants for reactant in self.reactants):
            self.energy_from_reactions = 0

        logger.debug("-== Phase 5. Distributing energy over reactants (with energy loss) ==-")
        energy_lost = self.dissipate_energy()
        self.system_energy -= energy_lost

        logger.debug("-== Phase 6. Removing strains  ==-")
        self.remove_strain()

        logger.info("*== Final reaction state ==*")
        logger.info(self.summary)
        logger.info(f"Energy lost: {energy_lost:8.2f} kJ/mol ({self.energy_loss*100:4.2f} %)")

        self.__update_reaction_summary(reactants=sorted(starting_reactants),
                                       products=sorted([x.symbol for x in self.reactants]))

        logger.info("*===============================*")
        return self.reactants, energy_lost


def __hydrogen_combustion_example():
    print("\n\n *===== Hydrogen Combustion =====*")
    h2_1 = Compound.from_formula("H2",0)
    h2_2 = Compound.from_formula("H2",0)
    o2 = Compound.from_formula("O2",0)

    h2_1.dissipate_energy(0.6)
    h2_2.dissipate_energy(0.6)
    o2.dissipate_energy(0.6)

    reaction = ReactionEngine(800)
    reaction.add_reactants([h2_1, h2_2, o2])
    products_first_reaction, _ = reaction.evaluate_reaction()
    products_second_reaction, _ = reaction.evaluate_reaction()

    return products_second_reaction

def __large_scale_hydrogen_combustion_example():
    comps = [Compound.from_formula("H2", 50) for _ in range(20)]
    comps.extend([Compound.from_formula("O2", 40) for _ in range(10)])

    reaction = ReactionEngine(0)
    reaction.add_reactants(comps)
    p1, _ = reaction.evaluate_reaction()
    p2, _ = reaction.evaluate_reaction()
    return p2

def __carbic_acid_synthesis_example():
    print("\n\n *===== Carbic Acid Synthesis =====*")
    co2 = Compound.from_formula("CO2", 0)
    h2o = Compound.from_formula("H2O", 0)

    reaction = ReactionEngine(0, energy_loss=0.2)
    reaction.add_reactant(co2)
    reaction.add_reactant(h2o)

    result, _ = reaction.evaluate_reaction()


if __name__ == "__main__":
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)

    # Example reactions handled by this module

    # __hydrogen_combustion_example()
    # __large_scale_hydrogen_combustion_example()
    __carbic_acid_synthesis_example()

    # def get_coh()->Compound:
    #     coh = Compound.from_formula("COH", 0)
    #     return coh
    #
    # coh = get_coh()
    #
    # reac = ReactionEngine()
    # reac.add_reactant(coh)
    # reac.add_reactant(Compound.atom('O'))
    # prod, _ = reac.evaluate_reaction()
    # prod[0].draw_compound()
    # prod[0].summary()
