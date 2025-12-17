from __future__ import annotations
from enum import Enum
from typing import List, Dict, Optional, Tuple


class Component:
    class FundamentalComponentNotDefined(Exception):
        pass

    fundamental_components : Dict[str, Component] = {}

    def __init__(self, mass: float, symbol: str, electrons_in_covalence: int, optimal_electrons: int):
        assert optimal_electrons >= electrons_in_covalence
        self.mass : float = mass
        self.symbol : str = symbol
        self.cov_electrons : int = electrons_in_covalence
        self.optimal_electrons : int = optimal_electrons

    def electrons_needed(self) -> int:
        return self.optimal_electrons - self.cov_electrons

    @staticmethod
    def get_fundamental(component_name: str) -> Component:
        try:
            return Component.fundamental_components[component_name.upper()]
        except KeyError:
            raise Component.FundamentalComponentNotDefined

    def __repr__(self):
        return f'{self.symbol:2s} Mass: {self.mass:3.2f} Electrons: {self.cov_electrons:1d}'


Component.fundamental_components = {
    'H': Component(mass=2, symbol="H", electrons_in_covalence=1, optimal_electrons=2),
    'C': Component(mass=12, symbol="C", electrons_in_covalence=4, optimal_electrons=8),
    'O': Component(mass=16, symbol="O", electrons_in_covalence=6, optimal_electrons=8),
}

class Bond:
    class BondNotDefinedError(Exception):
        pass

    '''Class defining energy values of bonds between two atoms

    Source: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Bond_Energies
    '''
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
        # Ensure alphabetical order
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
        self.bonds: List[Tuple[int, int, int]] = []

        mass = sum([x.mass for x in self.components])

        symbols = []
        symbols_count = {com[0]: len(com) for com in (list(x.symbol for x in self.components if x.symbol == y.symbol) for y in set(self.components))}
        list(map(lambda x: symbols.extend([str(x[0]), str(x[1])]), sorted(symbols_count.items(), key=lambda x: x[0], reverse=False)))
        symbol = "".join(symbols)
        optimal_electrons = sum([x.optimal_electrons for x in components])

        self.stable : bool = False
        electrons = sum([x.cov_electrons for x in self.components])
        if electrons // 8 == electrons / 8:
            self.stable = True

        super().__init__(mass, symbol, electrons, optimal_electrons)

        self.stable, self.remaining_energy = self.evaluate_bonds(provided_energy)
        print(self.remaining_energy)
        print(self.bonds)

    def evaluate_bonds(self, provided_energy: int):
        '''
        bonds = []

        available_components = self.components.copy()
        for m in range(0, 10):
            for bond in Bond.all_bonds():
                if len(available_components) < 2:
                    break
                try:
                    comp_A = next(filter(lambda x: x[1].symbol == bond.component_A.symbol, enumerate(available_components)))
                    comp_B = next(filter(lambda x: x[1].symbol == bond.component_B.symbol, enumerate(available_components)))
                except StopIteration:
                    continue
                if comp_A[0] == comp_B[0]:
                    continue


                to_remove = []
                for b in bonds:
                    if ((bond.component_A, bond.component_B) == (b.component_A, bond.component_B)
                            and ((provided_energy+b.energy-bond.energy)>= 0)
                            and (b.multiplicity < bond.multiplicity)):
                        provided_energy += b.energy
                        to_remove.append(b)
                for x in to_remove:
                    bonds.remove(x)

                if (provided_energy - bond.energy) < 0:
                    break
                provided_energy -= bond.energy

                available_components[comp_A[0]].cov_electrons += 1 * bond.multiplicity
                if available_components[comp_A[0]].cov_electrons >= 8:
                    available_components.pop(comp_A[0])
                    comp_B = (comp_B[0]-1, comp_B[1])

                available_components[comp_B[0]].cov_electrons += 1 * bond.multiplicity
                if available_components[comp_B[0]].cov_electrons >= 8:
                    available_components.pop(comp_B[0])

                bonds.append(Bond(comp_A[1], comp_B[1], bond.multiplicity))

        self.bonds = bonds

        for b in self.bonds:
            print(b)
        print(provided_energy)
        return provided_energy
        '''
        n = len(self.components)
        energy_used = 0

        max_iterations = 20

        for iteration in range(max_iterations):
            best_bond = None
            best_energy = float('inf')
            best_indices = None

            for i in range(n):
                for j in range(i + 1, n):
                    comp_i = self.components[i]
                    comp_j = self.components[j]

                    if comp_i.electrons_needed() == 0 and comp_j.electrons_needed() == 0:
                        continue

                    for mult in [1, 2, 3]:
                        bond_energy = Bond.get_bond_energy(comp_i.symbol, comp_j.symbol, mult)
                        if bond_energy is None:
                            continue

                        if energy_used + bond_energy > provided_energy:
                            continue

                        electrons_provided = mult
                        would_help = (comp_i.electrons_needed() >= electrons_provided or
                                      comp_j.electrons_needed() >= electrons_provided)

                        if not would_help:
                            continue

                        if bond_energy < best_energy:
                            best_energy = bond_energy
                            best_bond = mult
                            best_indices = (i, j)

            if best_bond is not None:
                i, j = best_indices
                mult = best_bond

                self.bonds.append((i, j, mult))
                self.components[i].cov_electrons += mult
                self.components[j].cov_electrons += mult
                energy_used += best_energy
            else:
                break

        all_stable = all(comp.cov_electrons >= comp.optimal_electrons for comp in self.components)

        return all_stable, provided_energy - energy_used


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
    comp = Compound([Component.get_fundamental('H'), Component.get_fundamental('H'), Component.get_fundamental('O')], 1200)
    print(comp)