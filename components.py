from __future__ import annotations

import itertools
from typing import List, Dict, Optional, Tuple, Generator

from scipy.special import binom
from tqdm import tqdm

class Component:
    class FundamentalComponentNotDefined(Exception):
        pass

    fundamental_components : Dict[str, Dict] = {
        'H': {'mass':2, 'symbol':'H', 'electrons_in_covalence':1, 'optimal_electrons':2},
        'C': {'mass':12, 'symbol':'C', 'electrons_in_covalence':4, 'optimal_electrons':8},
        'O': {'mass':16, 'symbol':'O', 'electrons_in_covalence':6, 'optimal_electrons':8},
    }

    def __init__(self, mass: float, symbol: str, electrons_in_covalence: int, optimal_electrons: int):
        assert optimal_electrons >= electrons_in_covalence
        self.mass : float = mass
        self.symbol : str = symbol
        self.cov_electrons : int = electrons_in_covalence
        self.optimal_electrons : int = optimal_electrons

    @property
    def electrons_needed(self) -> int:
        return self.optimal_electrons - self.cov_electrons

    @staticmethod
    def get_fundamental(component_name: str) -> Component:
        try:
            return Component(**Component.fundamental_components[component_name.upper()])
        except KeyError:
            raise Component.FundamentalComponentNotDefined

    def __repr__(self):
        return f'{self.symbol:2s} Mass: {self.mass:3.2f} Electrons: {self.cov_electrons:1d}'

def bond_combinations(components_bonds: Tuple[Tuple[int]]):
    number_of_bonds : int = int(binom(len(components_bonds),2))
    possible_bonds = list([0 for _ in range(number_of_bonds)])

    i = 0
    for a, comp_A in enumerate(components_bonds[:-1]):
        for b, comp_B in enumerate(components_bonds[a+1:]):
            possible_values = comp_A if len(comp_A) < len(comp_B) else comp_B
            possible_bonds[i] = possible_values
            i+=1

    return itertools.product(*possible_bonds)


class Bond:
    '''Class defining energy values of bonds between two atoms

    Source: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Bond_Energies
    '''

    class BondNotDefinedError(Exception):
        pass

    BondEnergy = {
        1: {
            'H-H': 432,
            'H-O': 467,
            'C-H': 413,
            'C-C': 347,
            'C-O': 358,
            'O-O': 146,
        },
        2: {
            'C-C': 614,
            'C-O': 745,
            'O-O': 495,
        },
        3: {
            'C-C': 839,
            'C-O': 1072,
        }
    }

    def __init__(self, component_A: Component, component_B: Component, multiplicity: int = 1):
        self.component_A : Component = component_A
        self.component_B : Component = component_B
        if self.component_A.symbol > self.component_B.symbol:
            temp = self.component_A
            self.component_A = self.component_B
            self.component_B = temp

        self.multiplicity : int = multiplicity
        try:
            self.energy = Bond.BondEnergy[multiplicity][f'{self.component_A.symbol}-{self.component_B.symbol}']
        except KeyError:
            raise Bond.BondNotDefinedError

    def __repr__(self):
        return f"{self.component_A.symbol:>2s}{['-','=','≡'][self.multiplicity-1]}{self.component_B.symbol:2s} Stored energy: {self.energy} kJ/mol"

    @staticmethod
    def get_bond_energy(symbol_A: str, symbol_B: str, multiplicity: int) -> Optional[int]:
        """Get bond energy if it exists, otherwise return None"""
        if symbol_A > symbol_B:
            symbol_A, symbol_B = symbol_B, symbol_A

        try:
            return Bond.BondEnergy[multiplicity][f'{symbol_A}-{symbol_B}']
        except KeyError:
            return None

    @staticmethod
    def all_bonds() -> List[Bond]:
        bonds = []
        for k, v in Bond.BondEnergy.items():
            for k2, v2 in v.items():
                bonds.append(Bond(*map(lambda x: Component.get_fundamental(x), k2.split('-')), multiplicity=k))
        return sorted(bonds, key=lambda x: x.energy, reverse=False)




class Compound(Component):
    def __init__(self, components: List[Component], provided_energy: int):
        self.components : List[Component] = sorted(components, key=lambda x: x.symbol, reverse=False)
        self.bonds : List[Tuple[Bond, int, int]] = []

        mass = sum([x.mass for x in self.components])

        symbols = []
        symbols_count = {com[0]: len(com) for com in (list(x.symbol for x in self.components if x.symbol == y.symbol) for y in set(self.components))}
        list(map(lambda x: symbols.extend([str(x[0]), str(x[1])]), sorted(symbols_count.items(), key=lambda x: x[0], reverse=False)))
        symbol = "".join(symbols)
        optimal_electrons = sum([x.optimal_electrons for x in components])

        self.stable : bool = False
        electrons = sum([x.cov_electrons for x in self.components])

        super().__init__(mass, symbol, electrons, optimal_electrons)

        self.stable, self.remaining_energy = (self.optimal_bonds(provided_energy))
        print(self.remaining_energy)
        print(self.bonds)

    def optimal_bonds(self, provided_energy: int):
        '''
        Function finding optimal bonds for compound given some energy for synthesising bonds

        Brute force method - evaluate all possible bonds combinations and choose one with best score
        (very inefficient : O(4^C(n,2)), where n - number of components)

        '''
        possible_bonds = tuple((tuple((i for i in range(min(3, max(1, c.electrons_needed))+1))) for c in self.components))
        possible_combinations = bond_combinations(possible_bonds)

        best_score = 0
        best_comb = None
        for comb in tqdm(possible_combinations):
            score = Compound.evaluate_bond(comb, tuple(self.components), provided_energy)
            if score > best_score:
                best_score = score
                best_comb = comb


        self.bonds, energy_used = Compound.combination_to_bonds(best_comb, self.components)

        all_stable = all(comp.cov_electrons >= comp.optimal_electrons for comp in self.components)

        return all_stable, provided_energy - energy_used

    @staticmethod
    def evaluate_bond(bonds: Tuple[int], components: Tuple[Component], max_energy: float) -> float:
        score : float = 0.0
        energy_used : int = 0
        a = 0
        b = 0
        bonded_components = [False for _ in range(len(components))]
        components_electrons = [x.cov_electrons for x in components]

        for i, bond in enumerate(bonds):
            b += 1
            if b >= len(components):
                a += 1
                b = a+1
            if bond == 0:
                continue
            comp_a = components[a]
            comp_b = components[b]
            components_electrons[a] += bond
            components_electrons[b] += bond
            if components_electrons[a] > components[a].optimal_electrons or \
                components_electrons[b] > components[b].optimal_electrons:
                return 0.0
            bonded_components[a] = True
            bonded_components[b] = True
            energy = Bond.get_bond_energy(comp_a.symbol, comp_b.symbol, int(bond))
            if energy is None:
                return 0.0
            energy_used += energy

        satisfied_components = [x == components[i].optimal_electrons for i,x in enumerate(components_electrons)]

        if energy_used > max_energy:
            return 0.0

        if all(bonded_components):
            score += max_energy
        if all(satisfied_components):
            score += max_energy

        score += max_energy - energy_used

        return score

    @staticmethod
    def combination_to_bonds(bond_combination: Tuple[int], components: List[Component]) -> Tuple[List[Tuple[Bond, int, int]], int]:
        bonds = []
        energy_used = 0
        a : int = 0
        b : int = 0
        for i, mult in enumerate(bond_combination):
            b += 1
            if b >= len(components):
                a += 1
                b = a + 1
            if mult == 0:
                continue
            components[a].cov_electrons += mult
            components[b].cov_electrons += mult
            energy_used += Bond.get_bond_energy(components[a].symbol, components[b].symbol, mult)
            bonds.append((Bond(components[a], components[b], mult), a, b))
        return bonds, energy_used

    def __repr__(self):
        return f'{self.symbol:12s} | STABLE: {self.stable:1d} | Mass: {self.mass:4.2f}'

    @staticmethod
    def synthesis(src: Component, other: Component, provided_energy: int) -> Compound:
        if isinstance(other, Compound):
            other = other.components
        else:
            other = [other]
        if isinstance(src, Compound):
            src = src.components
        else:
            src = [src]
        return Compound(src+other, provided_energy)


if __name__ == "__main__":
    # testing common compounds
    h2o = Compound([Component.get_fundamental('H'), Component.get_fundamental('O'), Component.get_fundamental('H'), Component.get_fundamental('O')], 20000)
    print(h2o)
    for el in h2o.components:
        print(el)
    assert h2o.stable == True

    ch4 = Compound([Component.get_fundamental('H'), Component.get_fundamental('C'), Component.get_fundamental('H'), Component.get_fundamental('H'), Component.get_fundamental('H')], 2000)
    print(ch4)
    for el in ch4.components:
        print(el)
    assert ch4.stable == True

    co2 = Compound([Component.get_fundamental('C'), Component.get_fundamental('O'), Component.get_fundamental('O')], 2000)
    print(co2)
    for el in co2.components:
        print(el)
    assert co2.stable == True

    comps = [Component.get_fundamental('H') for _ in range(5)]
    impossible = Compound(comps, 10000)
    print(impossible)
    for el in impossible.components:
        print(el)
    assert impossible.stable == False

    '''
    # complex compound 
    comps = [Component.get_fundamental('C') for _ in range(6)]
    comps.extend([Component.get_fundamental('O') for _ in range(6)])
    comps.extend([Component.get_fundamental('H') for _ in range(12)])
    glucose = Compound(comps, 500000)
    print(glucose)
    '''