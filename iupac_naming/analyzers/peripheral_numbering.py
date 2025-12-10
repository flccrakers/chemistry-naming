"""
Peripheral Numbering for Fused Ring Systems
============================================

Implements the peripheral numbering system for fused polycyclic systems
according to IUPAC 2013 Blue Book P-25.3.3.

The peripheral numbering is essential for:
1. Determining fusion locants
2. Identifying principal positions
3. Constructing systematic fusion names
"""

from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict, deque
from rdkit import Chem


class PeripheralNumberer:
    """
    Calculates peripheral numbering for fused ring systems.
    
    Reference: P-25.3.3
    
    The peripheral numbering follows these rules:
    1. Start at the uppermost rightmost atom
    2. Proceed clockwise around the periphery
    3. Interior atoms (shared by 3+ rings) are not numbered peripherally
    4. Fusion atoms get sequential numbers
    """
    
    def __init__(self, mol: Chem.Mol, fused_atoms: Set[int]):
        """
        Initialize the numberer.
        
        Args:
            mol: RDKit molecule
            fused_atoms: Set of atom indices in the fused system
        """
        self.mol = mol
        self.fused_atoms = fused_atoms
        self.coords = {}  # atom_idx -> (x, y)
        self._extract_coordinates()
        
        # Results
        self.peripheral_order: List[int] = []
        self.peripheral_numbers: Dict[int, int] = {}  # atom_idx -> peripheral number
        self.interior_atoms: Set[int] = set()
        self.fusion_bonds: List[Tuple[int, int]] = []
    
    def _extract_coordinates(self):
        """Extract 2D coordinates from molecule."""
        conf = self.mol.GetConformer() if self.mol.GetNumConformers() > 0 else None
        
        for atom_idx in self.fused_atoms:
            if conf:
                pos = conf.GetAtomPosition(atom_idx)
                self.coords[atom_idx] = (pos.x, pos.y)
            else:
                # No coordinates - generate approximate ones
                self.coords[atom_idx] = (0, 0)
    
    def compute_numbering(self) -> Dict[int, int]:
        """
        Compute the peripheral numbering.
        
        Returns:
            Dictionary mapping atom index to peripheral number
        """
        # Step 1: Find peripheral vs interior atoms
        self._classify_atoms()
        
        if not self.peripheral_order:
            return {}
        
        # Step 2: Find starting atom (uppermost rightmost)
        start_atom = self._find_start_atom()
        
        # Step 3: Trace periphery clockwise
        self._trace_periphery(start_atom)
        
        # Step 4: Assign numbers
        for i, atom_idx in enumerate(self.peripheral_order):
            self.peripheral_numbers[atom_idx] = i + 1
        
        return self.peripheral_numbers
    
    def _classify_atoms(self):
        """
        Classify atoms as peripheral or interior.
        
        Interior atoms are those that are shared by 3+ rings or
        are completely surrounded by other fused atoms.
        """
        # Count how many rings each atom belongs to
        ring_info = self.mol.GetRingInfo()
        atom_ring_count = defaultdict(int)
        
        for ring in ring_info.AtomRings():
            ring_set = set(ring)
            if ring_set & self.fused_atoms:
                for atom in ring:
                    if atom in self.fused_atoms:
                        atom_ring_count[atom] += 1
        
        # Atoms in 3+ rings are interior
        for atom, count in atom_ring_count.items():
            if count >= 3:
                self.interior_atoms.add(atom)
        
        # Also check connectivity - if all neighbors are in fused system,
        # and atom has no external bonds, it might be interior
        for atom_idx in self.fused_atoms:
            atom = self.mol.GetAtomWithIdx(atom_idx)
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            
            # Count neighbors in vs out of fused system
            internal_neighbors = sum(1 for n in neighbors if n in self.fused_atoms)
            external_neighbors = sum(1 for n in neighbors if n not in self.fused_atoms)
            
            # Atom with no external bonds and high internal connectivity might be interior
            if external_neighbors == 0 and internal_neighbors >= 3:
                # But only if it's in multiple rings
                if atom_ring_count[atom_idx] >= 2:
                    self.interior_atoms.add(atom_idx)
        
        # Peripheral atoms are those not interior
        peripheral = self.fused_atoms - self.interior_atoms
        self.peripheral_order = list(peripheral)
    
    def _find_start_atom(self) -> int:
        """
        Find the starting atom for numbering.
        
        Reference: P-25.3.3.1
        
        The starting atom is:
        1. The uppermost atom
        2. Among those, the rightmost
        """
        if not self.peripheral_order:
            return -1
        
        peripheral = self.fused_atoms - self.interior_atoms
        
        # Sort by y (descending = uppermost) then x (descending = rightmost)
        def sort_key(atom_idx):
            x, y = self.coords.get(atom_idx, (0, 0))
            return (-y, -x)  # Negative for descending order
        
        sorted_atoms = sorted(peripheral, key=sort_key)
        return sorted_atoms[0] if sorted_atoms else -1
    
    def _trace_periphery(self, start: int):
        """
        Trace around the periphery clockwise from start atom.
        
        Uses the righthand rule: at each atom, choose the bond
        that makes the smallest clockwise angle from the incoming direction.
        """
        if start < 0:
            return
        
        peripheral = self.fused_atoms - self.interior_atoms
        
        # Build adjacency for peripheral atoms within fused system
        adj = defaultdict(list)
        for atom_idx in peripheral:
            atom = self.mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx in peripheral:
                    adj[atom_idx].append(n_idx)
        
        # Trace using angle-based selection
        visited = set()
        order = []
        current = start
        prev_direction = (0, -1)  # Initially coming from "above" (start going down)
        
        max_iterations = len(peripheral) * 2  # Safety limit
        iterations = 0
        
        while iterations < max_iterations:
            iterations += 1
            
            if current in visited:
                if current == start and len(order) > 2:
                    break  # Complete loop
                # Try to continue if possible
                unvisited = [n for n in adj[current] if n not in visited]
                if not unvisited:
                    break
            
            visited.add(current)
            order.append(current)
            
            # Get neighbors
            neighbors = [n for n in adj[current] if n != order[-2] if len(order) > 1]
            if len(order) == 1:
                neighbors = adj[current]
            
            if not neighbors:
                break
            
            # Select next atom: the one making smallest clockwise angle
            current_pos = self.coords.get(current, (0, 0))
            next_atom = self._select_clockwise(current, neighbors, prev_direction)
            
            if next_atom is None:
                break
            
            # Update direction
            next_pos = self.coords.get(next_atom, (0, 0))
            prev_direction = (next_pos[0] - current_pos[0], 
                            next_pos[1] - current_pos[1])
            
            current = next_atom
        
        self.peripheral_order = order
    
    def _select_clockwise(self, current: int, neighbors: List[int], 
                          incoming_dir: Tuple[float, float]) -> Optional[int]:
        """
        Select the neighbor that makes the smallest clockwise angle.
        
        Args:
            current: Current atom index
            neighbors: List of neighbor atom indices
            incoming_dir: Direction vector we came from
        
        Returns:
            Index of best neighbor, or None
        """
        import math
        
        if not neighbors:
            return None
        
        if len(neighbors) == 1:
            return neighbors[0]
        
        current_pos = self.coords.get(current, (0, 0))
        
        # Normalize incoming direction and get reference angle
        in_len = math.sqrt(incoming_dir[0]**2 + incoming_dir[1]**2)
        if in_len == 0:
            in_len = 1
        
        ref_angle = math.atan2(incoming_dir[1], incoming_dir[0])
        
        best_neighbor = None
        best_angle = float('inf')
        
        for n in neighbors:
            n_pos = self.coords.get(n, (0, 0))
            dx = n_pos[0] - current_pos[0]
            dy = n_pos[1] - current_pos[1]
            
            out_angle = math.atan2(dy, dx)
            
            # Calculate clockwise angle from incoming direction
            # Clockwise in screen coords (y-down) = counterclockwise in math coords
            angle_diff = ref_angle - out_angle
            
            # Normalize to [0, 2π)
            while angle_diff < 0:
                angle_diff += 2 * math.pi
            while angle_diff >= 2 * math.pi:
                angle_diff -= 2 * math.pi
            
            # We want the smallest positive clockwise angle (but not 0)
            if angle_diff < 0.01:  # Nearly same direction (continuing straight)
                angle_diff = 2 * math.pi  # Put it last
            
            if angle_diff < best_angle:
                best_angle = angle_diff
                best_neighbor = n
        
        return best_neighbor
    
    def get_face_letters(self, component_atoms: Set[int]) -> Dict[Tuple[int, int], str]:
        """
        Assign face letters to fusion bonds.
        
        Reference: P-25.3.2.4
        
        Faces are labeled with lowercase letters (a, b, c...)
        starting from position 1-2 and going around.
        
        Args:
            component_atoms: Atoms of a fused component
        
        Returns:
            Dictionary mapping (atom1, atom2) to face letter
        """
        if not self.peripheral_numbers:
            self.compute_numbering()
        
        face_letters = {}
        letter_ord = ord('a')
        
        # Find bonds between consecutive peripheral positions
        for i in range(len(self.peripheral_order)):
            a1 = self.peripheral_order[i]
            a2 = self.peripheral_order[(i + 1) % len(self.peripheral_order)]
            
            # Check if there's a bond
            bond = self.mol.GetBondBetweenAtoms(a1, a2)
            if bond:
                num1 = self.peripheral_numbers.get(a1, 0)
                num2 = self.peripheral_numbers.get(a2, 0)
                
                # Bond between positions n and n+1 is face letter
                letter = chr(letter_ord)
                face_letters[(a1, a2)] = letter
                face_letters[(a2, a1)] = letter  # Both directions
                letter_ord += 1
        
        return face_letters
    
    def get_fusion_locants(self, attached_atoms: Set[int], 
                          base_atoms: Set[int]) -> Tuple[List[int], List[int]]:
        """
        Get fusion locants between an attached component and base.
        
        Reference: P-25.3.2
        
        Returns:
            (attached_locants, base_locants)
        """
        if not self.peripheral_numbers:
            self.compute_numbering()
        
        # Find shared atoms
        shared = attached_atoms & base_atoms
        
        if len(shared) < 2:
            return ([], [])
        
        # Get peripheral numbers for shared atoms
        shared_with_numbers = []
        for atom in shared:
            if atom in self.peripheral_numbers:
                shared_with_numbers.append((atom, self.peripheral_numbers[atom]))
        
        # Sort by peripheral number
        shared_with_numbers.sort(key=lambda x: x[1])
        
        # Return locant lists
        locants = [num for atom, num in shared_with_numbers]
        
        return (locants, locants)  # Same for both if using base numbering


def compute_fusion_descriptor(mol: Chem.Mol,
                             attached_component: Dict,
                             base_component: Dict,
                             all_fused_atoms: Set[int],
                             prime_level: int = 0) -> str:
    """
    Compute the fusion descriptor string for an attached component.
    
    Reference: P-25.3.2
    
    Args:
        mol: RDKit molecule
        attached_component: Dict with 'name', 'atoms', etc.
        base_component: Dict with base component info
        all_fused_atoms: All atoms in the fused system
        prime_level: Number of primes for locants (0, 1, 2...)
    
    Returns:
        Fusion descriptor like "[4',5':5,6]" or "[2,3-b]"
    """
    # Initialize numberer
    numberer = PeripheralNumberer(mol, all_fused_atoms)
    numberer.compute_numbering()
    
    # Find shared atoms
    shared = attached_component['atoms'] & base_component['atoms']
    
    if len(shared) < 2:
        return ""
    
    # Get locants
    attached_locants, base_locants = numberer.get_fusion_locants(
        attached_component['atoms'], base_component['atoms']
    )
    
    if not attached_locants:
        return ""
    
    # Format with primes
    def format_locant(num: int, primes: int) -> str:
        return str(num) + "'" * primes
    
    attached_str = ','.join(format_locant(loc, prime_level) for loc in attached_locants)
    base_str = ','.join(str(loc) for loc in base_locants)
    
    return f"[{attached_str}:{base_str}]"


def format_fusion_name(components: List[Dict], 
                       base: Dict,
                       mol: Chem.Mol,
                       bridges: List[Dict] = None) -> str:
    """
    Format a complete fusion name with all locants.
    
    Reference: P-25.3.4
    
    Args:
        components: List of attached component dicts
        base: Base component dict
        mol: RDKit molecule
        bridges: Optional list of bridge dicts
    
    Returns:
        Complete fusion name like "pyrido[1',2':1,2]imidazo[4,5:5,6]pyrazino[2,3-b]phenazine"
    """
    # Get all fused atoms
    all_atoms = set()
    all_atoms.update(base['atoms'])
    for comp in components:
        all_atoms.update(comp['atoms'])
    
    # Order components by distance from base (outermost first)
    ordered = order_components_by_distance(components, base)
    
    # Build name parts
    parts = []
    
    # Add bridges first
    if bridges:
        for bridge in bridges:
            bridge_name = bridge.get('bridge_type', 'unknown')
            loc1 = bridge.get('locant1', '?')
            loc2 = bridge.get('locant2', '?')
            parts.append(f"{loc1},{loc2}-{bridge_name}")
    
    # Add fusion components
    prime_level = len(ordered) - 1  # Outermost gets most primes
    
    for comp in ordered:
        prefix = comp.get('fusion_prefix', comp['name'] + 'o')
        
        # Calculate fusion descriptor
        descriptor = compute_fusion_descriptor(
            mol, comp, base, all_atoms, prime_level
        )
        
        parts.append(f"{prefix}{descriptor}")
        prime_level -= 1
    
    # Add base name
    parts.append(base['name'])
    
    return ''.join(parts)


def order_components_by_distance(components: List[Dict], base: Dict) -> List[Dict]:
    """
    Order attached components from outermost to innermost.
    
    Reference: P-25.3.4.2.3
    
    Components directly attached to base come last.
    Components attached to those come before, etc.
    """
    if not components:
        return []
    
    def shares_atoms(c1, c2):
        return len(c1['atoms'] & c2['atoms']) > 0
    
    # Separate by direct attachment to base
    direct = [c for c in components if shares_atoms(c, base)]
    indirect = [c for c in components if not shares_atoms(c, base)]
    
    # BFS to order indirect
    ordered_indirect = []
    remaining = indirect.copy()
    current_level = direct
    
    while remaining:
        next_level = []
        for comp in remaining:
            for attached in current_level:
                if shares_atoms(comp, attached):
                    next_level.append(comp)
                    break
        
        if next_level:
            # Prepend to get outermost first
            ordered_indirect = next_level + ordered_indirect
            for c in next_level:
                if c in remaining:
                    remaining.remove(c)
            current_level = next_level
        else:
            break
    
    # Final order: outermost indirect, then direct
    return ordered_indirect + direct
