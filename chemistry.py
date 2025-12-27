from __future__ import annotations
from typing import List, Dict, Optional, Tuple, Generator
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


class Atom:
    class AtomNotDefined(Exception):
        pass

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

    @property
    def electrons_needed(self) -> int:
        return self.optimal_electrons - self.cov_electrons

    @staticmethod
    def get(atom_symbol: str) -> Atom:
        try:
            return Atom(**Atom.basic_atoms[atom_symbol.upper()])
        except KeyError:
            raise Atom.AtomNotDefined

    def __repr__(self):
        return f'{self.symbol:2s} Mass: {self.mass:3.2f} Electrons: {self.cov_electrons:1d}'


class Bond:
    '''
    Class defining energy values of bonds between two atoms
    Source: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Bond_Energies
    '''
    class BondNotDefinedError(Exception):
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

    def __repr__(self):
        return f"{self.component_A.symbol:>2s}{['-', '=', '≡'][self.multiplicity - 1]}{self.component_B.symbol:2s} Stored energy: {self.energy} kJ/mol"

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
        self.components = components
        self.bonds = []
        self.stable = False
        self.preserve_bonds = preserve_bonds or []
        self.remaining_energy = provided_energy

        self.mass = sum([x.mass for x in self.components])

        symbols = []
        symbols_count = {com[0]: len(com) for com in (list(x.symbol for x in self.components if x.symbol == y.symbol) for y in set(self.components))}
        list(map(lambda x: symbols.extend([str(x[0]), str(x[1])]), sorted(symbols_count.items(), key=lambda x: x[0], reverse=False)))
        self.symbol = "".join(symbols)

        self.optimal_electrons = sum([x.optimal_electrons for x in components])
        self.electrons = sum([x.cov_electrons for x in self.components])

    @staticmethod
    def atom(atom_symbol:str, provided_energy:int = 0) -> Compound:
        return Compound([Atom.get(atom_symbol)], provided_energy=provided_energy)

    def __repr__(self):
        return f'{self.symbol:12s} | STABLE: {self.stable:1d} | Mass: {self.mass:4.2f}'

    def show_structure(self):
        print(f"\n{self.symbol} Structure:")
        print(f"Atoms ({len(self.components)}):")
        for i, comp in enumerate(self.components):
            status = "✓" if comp.cov_electrons >= comp.optimal_electrons else "✗"
            print(f"  [{i}] {comp.symbol}: {comp.cov_electrons}/{comp.optimal_electrons} electrons {status}")

        print(f"\nBonds ({len(self.bonds)}):")
        for bond, i, j in self.bonds:
            bond_symbol = ['-', '=', '≡'][bond.multiplicity - 1]
            print(
                f"  [{i}]{self.components[i].symbol} {bond_symbol} {self.components[j].symbol}[{j}] ({bond.energy} kJ/mol)")

        print(f"\nStable: {self.stable}, Energy remaining: {self.remaining_energy} kJ/mol")

    def draw_compound(self):
        import networkx as nx
        import matplotlib.pyplot as plt
        print(self.bonds)

        single = list([[f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$'] for bond, i, j in self.bonds if bond.multiplicity >= 1])
        double = list([[f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$'] for bond, i, j in self.bonds if bond.multiplicity >= 2])
        triple = list([[f'${self.components[i].symbol}_'+'{'f'{i}'+'}$', f'${self.components[j].symbol}_'+'{'f'{j}'+'}$'] for bond, i, j in self.bonds if bond.multiplicity >= 3])

        graph = nx.Graph()
        graph.add_edges_from(single)
        graph.add_edges_from(double)
        graph.add_edges_from(triple)

        pos = nx.spring_layout(graph, seed=42)

        nx.draw_networkx_nodes(graph, pos, node_size=600, node_color='white', edgecolors="black")
        nx.draw_networkx_labels(graph, pos, font_size=10)
        nx.draw_networkx_edges(graph, pos, arrows=True,edgelist=single, connectionstyle='arc3, rad = 0.1')
        nx.draw_networkx_edges(graph, pos,arrows=True, edgelist=double, connectionstyle='arc3, rad = 0.25')
        nx.draw_networkx_edges(graph, pos,arrows=True, edgelist=triple, connectionstyle='arc3, rad = 0.4')

        plt.show()

    def break_bonds(self, bond_indices: List[int]) -> int:
        """
        Break specified bonds and return energy released
        bond_indices: list of indices in self.bonds to break
        Returns: energy released (positive value)
        """
        energy_released = 0
        bonds_to_remove = []

        for idx in sorted(bond_indices, reverse=True):
            if idx < len(self.bonds):
                bond, i, j = self.bonds[idx]

                # Remove electrons from components
                self.components[i].cov_electrons -= bond.multiplicity
                self.components[j].cov_electrons -= bond.multiplicity

                # Energy is released when breaking bonds
                energy_released += bond.energy
                bonds_to_remove.append(idx)

        # Remove bonds
        for idx in bonds_to_remove:
            del self.bonds[idx]

        # Check if still stable
        self.stable = all(comp.cov_electrons >= comp.optimal_electrons
                          for comp in self.components)

        return energy_released

    @staticmethod
    def synthesize(comp_a: Compound,
                   comp_b: Compound,
                   provided_energy: int,
                   break_bonds_a: Optional[List[int]] = None,
                   break_bonds_b: Optional[List[int]] = None) -> Compound:
        """
        Synthesize a new compound from two reactants

        Args:
            comp_a: First reactant compound
            comp_b: Second reactant compound
            provided_energy: Energy provided for reaction
            break_bonds_a: Indices of bonds to break in reactant_a (optional)
            break_bonds_b: Indices of bonds to break in reactant_b (optional)

        Returns:
            New compound with reconfigured bonds
        """
        # Deep copy reactants to avoid modifying originals
        import copy
        r_a = copy.deepcopy(comp_a)
        r_b = copy.deepcopy(comp_b)

        # Break specified bonds and collect energy
        energy_from_breaking = 0
        if break_bonds_a:
            energy_from_breaking += r_a.break_bonds(break_bonds_a)
        if break_bonds_b:
            energy_from_breaking += r_b.break_bonds(break_bonds_b)

        # Total energy available
        total_energy = provided_energy + energy_from_breaking + comp_b.remaining_energy + comp_a.remaining_energy

        logger.info(f"=== Synthesis Reaction ===")
        logger.info(f"Breaking bonds released: {energy_from_breaking} kJ/mol")
        logger.info(f"Energy provided: {provided_energy} kJ/mol")
        logger.info(f"Total energy available: {total_energy} kJ/mol")

        # Combine components from both reactants
        combined_components = []
        atom_offset_b = len(r_a.components)

        # Track which bonds to preserve
        preserved_bonds = []

        # Add components from reactant A with current electron state
        for idx, comp in enumerate(r_a.components):
            new_comp = Atom(
                comp.mass, comp.symbol, comp.cov_electrons, comp.optimal_electrons
            )
            combined_components.append(new_comp)

        # Preserve ALL remaining bonds from reactant A (internal structure)
        for bond, i, j in r_a.bonds:
            preserved_bonds.append((i, j, bond.multiplicity))

        # Add components from reactant B with current electron state
        for idx, comp in enumerate(r_b.components):
            new_comp = Atom(
                comp.mass, comp.symbol, comp.cov_electrons, comp.optimal_electrons
            )
            combined_components.append(new_comp)

        # Preserve ALL remaining bonds from reactant B (internal structure)
        for bond, i, j in r_b.bonds:
            new_i = i + atom_offset_b
            new_j = j + atom_offset_b
            preserved_bonds.append((new_i, new_j, bond.multiplicity))

        logger.info(f"Preserved {len(preserved_bonds)} bonds from reactants")
        logger.info(f"Reactant A: atoms 0-{atom_offset_b - 1}")
        logger.info(f"Reactant B: atoms {atom_offset_b}-{len(combined_components) - 1}")

        product = Compound.__create_from_parts(
            combined_components,
            preserved_bonds,
            total_energy,
            atom_offset_b  # Boundary between reactant A and B
        )

        return product

    @staticmethod
    def __create_from_parts(components: List[Atom],
                            preserved_bonds: List[Tuple[int, int, int]],
                            provided_energy: int,
                            boundary_idx: int) -> 'Compound':
        """
        Create compound with directed bonding between two reactant groups
        Only forms bonds BETWEEN groups, not within them
        """
        product = Compound(components=components, provided_energy=provided_energy, preserve_bonds=preserved_bonds)

        # Calculate compound properties
        product.mass = sum(x.mass for x in components)
        symbols_count = {}
        for comp in components:
            symbols_count[comp.symbol] = symbols_count.get(comp.symbol, 0) + 1

        symbol_parts = []
        for sym in sorted(symbols_count.keys()):
            symbol_parts.extend([sym, str(symbols_count[sym])])
        product.symbol = "".join(symbol_parts)

        product.optimal_electrons = sum(x.optimal_electrons for x in components)

        # Perform directed synthesis
        product.stable, product.remaining_energy = product.__directed_bond_formation(boundary_idx)

        logger.info(f"Remaining energy: {product.remaining_energy}")
        logger.info(f"Bonds: {product.bonds}")

        return product

    def __directed_bond_formation(self, boundary_idx: int) -> Tuple[bool, int]:
        """
        Form bonds only between two reactant groups
        Preserves internal structure of each reactant
        PRIORITIZES creating a connected molecular structure

        Args:
            boundary_idx: Index separating reactant A (0 to boundary-1) from B (boundary to n-1)
        """
        n = len(self.components)
        electrons = [comp.cov_electrons for comp in self.components]
        optimal_electrons = [comp.optimal_electrons for comp in self.components]

        bonds_made = []
        energy_used = 0

        # First, add all preserved bonds (internal structures)
        preserved_set = set()
        for i, j, mult in self.preserve_bonds:
            if i < n and j < n:
                preserved_set.add((min(i, j), max(i, j)))
                bonds_made.append((i, j, mult, 0))  # 0 energy cost

        logger.info(f"Looking for reactive sites:")
        logger.info(f"Reactant A atoms (0-{boundary_idx - 1}):")
        for i in range(boundary_idx):
            need = optimal_electrons[i] - electrons[i]
            if need > 0:
                logger.info(f"  [{i}] {self.components[i].symbol}: needs {need} electrons")

        logger.info(f"Reactant B atoms ({boundary_idx}-{n - 1}):")
        for i in range(boundary_idx, n):
            need = optimal_electrons[i] - electrons[i]
            if need > 0:
                logger.info(f"  [{i}] {self.components[i].symbol}: needs {need} electrons")

        # Build connectivity graph from preserved bonds
        def get_connected_components():
            """Returns list of sets, each set contains indices of connected atoms"""
            parent = list(range(n))

            def find(x):
                if parent[x] != x:
                    parent[x] = find(parent[x])
                return parent[x]

            def union(x, y):
                px, py = find(x), find(y)
                if px != py:
                    parent[px] = py

            # Build union-find from existing bonds
            for i, j, mult, energy in bonds_made:
                union(i, j)

            # Group atoms by connected component
            components = {}
            for i in range(n):
                root = find(i)
                if root not in components:
                    components[root] = set()
                components[root].add(i)

            return list(components.values())

        # Check if structure is connected
        def is_connected():
            components = get_connected_components()
            return len(components) == 1

        logger.info(f"PHASE 1: Ensuring connectivity between reactants")

        connected_components = get_connected_components()
        logger.info(f"Initial connected components: {len(connected_components)}")

        # Find which components contain reactant A and B atoms
        component_a = None
        component_b = None
        for comp_set in connected_components:
            if any(i < boundary_idx for i in comp_set):
                component_a = comp_set
            if any(i >= boundary_idx for i in comp_set):
                component_b = comp_set

        # If reactants are in separate components, we MUST connect them
        if component_a != component_b:
            logger.info("Reactants are disconnected - finding bridge bond...")

            # Find the best bridge bond to connect the two reactants
            best_bridge = None
            best_bridge_score = -float('inf')

            for i in range(boundary_idx):  # Atoms from reactant A
                if i not in component_a:
                    continue

                for j in range(boundary_idx, n):  # Atoms from reactant B
                    if j not in component_b:
                        continue

                    need_i = optimal_electrons[i] - electrons[i]
                    need_j = optimal_electrons[j] - electrons[j]

                    if need_i <= 0 or need_j <= 0:
                        continue

                    # Try different bond multiplicities
                    for mult in range(1, min(need_i, need_j, 3) + 1):
                        energy = Bond.get_bond_energy(
                            self.components[i].symbol,
                            self.components[j].symbol,
                            mult
                        )

                        if energy and energy_used + energy <= self.remaining_energy:
                            # Score: prioritize satisfying electron needs
                            satisfaction = mult / need_i + mult / need_j
                            score = satisfaction * 1000 - energy

                            if score > best_bridge_score:
                                best_bridge_score = score
                                best_bridge = (i, j, mult, energy)

            # Add the bridge bond
            if best_bridge:
                i, j, mult, energy = best_bridge
                bonds_made.append((i, j, mult, energy))
                electrons[i] += mult
                electrons[j] += mult
                energy_used += energy
                preserved_set.add((min(i, j), max(i, j)))

                bond_symbol = ['-', '=', '≡'][mult - 1]
                logger.info(
                    f"  BRIDGE BOND: [{i}]{self.components[i].symbol} {bond_symbol} {self.components[j].symbol}[{j}] ({energy} kJ/mol)")
            else:
                logger.warning("No valid bridge bond found!")
        else:
            logger.info("Reactants already connected through preserved bonds")

        logger.info(f"PHASE 2: Satisfying remaining electron needs")

        # Find all possible bonds (excluding already bonded pairs)
        bonds = []

        for i in range(boundary_idx):  # Atoms from reactant A
            for j in range(boundary_idx, n):  # Atoms from reactant B
                if (min(i, j), max(i, j)) in preserved_set:
                    continue  # Already bonded

                need_i = optimal_electrons[i] - electrons[i]
                need_j = optimal_electrons[j] - electrons[j]

                if need_i <= 0 or need_j <= 0:
                    continue

                # Try different bond multiplicities
                for mult in range(1, min(need_i, need_j, 3) + 1):
                    energy = Bond.get_bond_energy(
                        self.components[i].symbol,
                        self.components[j].symbol,
                        mult
                    )

                    if energy and energy_used + energy <= self.remaining_energy:
                        # Calculate priority: prefer satisfying more needs
                        satisfaction = mult / need_i + mult / need_j
                        bonds.append((i, j, mult, energy, satisfaction))

        # Sort by satisfaction (higher is better), then by energy (lower is better)
        bonds.sort(key=lambda x: (-x[4], x[3]))

        logger.info(f"Found {len(bonds)} additional possible bonds")

        # Greedily add bonds
        bonds_added = 0
        for i, j, mult, energy, satisfaction in bonds:
            # Check if we can still add this bond
            if (electrons[i] + mult <= optimal_electrons[i] and
                    electrons[j] + mult <= optimal_electrons[j] and
                    energy_used + energy <= self.remaining_energy):
                bonds_made.append((i, j, mult, energy))
                electrons[i] += mult
                electrons[j] += mult
                energy_used += energy
                preserved_set.add((min(i, j), max(i, j)))
                bonds_added += 1

                bond_symbol = ['-', '=', '≡'][mult - 1]
                logger.info(
                    f"  Forming: [{i}]{self.components[i].symbol} {bond_symbol} {self.components[j].symbol}[{j}] ({energy} kJ/mol)")

        # Final connectivity check
        if is_connected():
            logger.info(f"Product is fully connected")
        else:
            connected_components = get_connected_components()
            logger.warning(f"Product has {len(connected_components)} disconnected fragments")
            for idx, comp_set in enumerate(connected_components):
                atoms = [f"{self.components[i].symbol}[{i}]" for i in sorted(comp_set)]
                logger.warning(f"  Fragment {idx + 1}: {' '.join(atoms)}")

        # Convert to Bond objects
        if bonds_made:
            self.bonds = []
            for i, j, mult, energy in bonds_made:
                bond_obj = Bond(self.components[i], self.components[j], mult)
                self.bonds.append((bond_obj, i, j))

            # Update component electron counts
            for k in range(n):
                self.components[k].cov_electrons = electrons[k]

            all_satisfied = all(
                electrons[k] >= optimal_electrons[k]
                for k in range(n)
            )

            return all_satisfied, self.remaining_energy - energy_used

        return False, self.remaining_energy

    @staticmethod
    def from_formula(formula: str, provided_energy: int) -> 'Compound':
        """
        Create a compound from a chemical formula string
        Example: from_formula("H2O", 2000) or from_formula("C6H12O6", 50000)
        """
        import re

        # Parse formula: C6H12O6 -> [('C', 6), ('H', 12), ('O', 6)]
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)

        compound = None
        for symbol, count in matches:
            count = int(count) if count else 1
            for _ in range(count):
                if compound is None:
                    compound = Compound.atom(symbol, provided_energy)
                else:
                    compound = Compound.synthesize(compound, Compound.atom(symbol), 0)
        return compound


def __glucose_synthesis():
    # Example synthesis of glucose molecule C6H12O6
    def get_coh()->Compound:
        oh = Compound.from_formula("OH", 2000)
        coh = Compound.synthesize(oh, Compound.atom("C"), 200)
        return coh

    glucose = Compound.from_formula("CHO", 2000)
    for i in range(5):
        glucose = Compound.synthesize(glucose, get_coh(), 2000)
        glucose = Compound.synthesize(glucose, Compound.atom("H"), 2000)
    glucose = Compound.synthesize(glucose, Compound.atom("H"), 2000)

    return glucose

def __ethane_synthesis() -> Compound:
    # Example synthesis of ethane C2H6
    def ch3() -> Compound:
        return Compound.from_formula("CH3", 2000)
    return Compound.synthesize(ch3(), ch3(), 100)


if __name__ == "__main__":
    logger.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)

    ethane = __ethane_synthesis()
    ethane.show_structure()
    ethane.draw_compound()

    glucose = __glucose_synthesis()
    glucose.show_structure()
    glucose.draw_compound()



