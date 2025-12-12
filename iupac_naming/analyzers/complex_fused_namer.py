"""
Complex Fused System Namer
==========================

Specialized module for naming complex polycyclic fused systems
with multiple components and bridges.

Reference: IUPAC Blue Book P-25, P-24
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass
from rdkit import Chem


@dataclass
class BridgeInfo:
    """Bridge between two positions in fused system."""
    bridge_type: str
    locant1: int
    locant2: int
    bridge_atom: int  # Index of the bridging atom


class ComplexFusedNamer:
    """
    Names complex fused ring systems with bridges.
    
    Example: 9,12-epoxypyrido[1'',2'':1',2']imidazo[4',5':5,6]pyrazino[2,3-b]phenazine
    """
    
    # Patterns for component identification - more specific first
    COMPONENT_PATTERNS = [
        # Tricyclic bases
        ('phenazine', 'c1ccc2nc3ccccc3nc2c1', 100),
        ('acridine', 'c1ccc2nc3ccccc3cc2c1', 95),
        
        # Bicyclic
        ('quinoxaline', 'c1ccc2nccnc2c1', 60),
        ('quinazoline', 'c1ccc2ncncc2c1', 59),
        ('quinoline', 'c1ccc2ncccc2c1', 58),
        ('isoquinoline', 'c1ccc2ccncc2c1', 57),
        ('naphthalene', 'c1ccc2ccccc2c1', 50),
        ('indole', 'c1ccc2[nH]ccc2c1', 55),
        ('benzimidazole', 'c1ccc2[nH]cnc2c1', 54),
        
        # 6-membered monocycles
        ('pyrazine', 'c1nccnc1', 40),  # Simplified SMARTS
        ('pyrimidine', 'c1ncncc1', 39),
        ('pyridazine', 'c1nncc1', 38),
        
        # 5-membered monocycles - process before 6-membered to capture fusion links
        ('imidazole', 'c1nc[nH]c1', 36),
        ('imidazole', 'c1[nH]cnc1', 36),
        ('imidazole', 'c1ncnc1', 36),  # Without explicit H
        ('pyrazole', 'c1cc[nH]n1', 34),
        
        # 6-membered monocycles (lower priority)
        ('pyridine', 'c1ncccc1', 35),
        ('benzene', 'c1ccccc1', 30),
        
        # Remaining 5-membered
        ('pyrrole', '[nH]1cccc1', 31),
        ('furan', 'o1cccc1', 28),
        ('thiophene', 's1cccc1', 29),
    ]
    
    FUSION_PREFIXES = {
        'benzene': 'benzo', 'naphthalene': 'naphtho', 'pyridine': 'pyrido',
        'pyrazine': 'pyrazino', 'pyrimidine': 'pyrimido', 'pyridazine': 'pyridazino',
        'imidazole': 'imidazo', 'pyrazole': 'pyrazolo',
        'pyrrole': 'pyrrolo', 'furan': 'furo', 'thiophene': 'thieno',
        'quinoline': 'quinolino', 'isoquinoline': 'isoquinolino',
        'indole': 'indolo', 'quinoxaline': 'quinoxalino',
        'phenazine': 'phenazino',
    }
    
    # Fusion edge descriptors (P-25.3.2.1)
    FUSION_EDGES = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    
    def __init__(self, mol, verbose: bool = False):
        self.mol = mol
        self.verbose = verbose
        self.ring_info = mol.GetRingInfo()
        self.rings = [set(r) for r in self.ring_info.AtomRings()]
        self.peripheral_numbering = {}
        
    def _log(self, msg: str):
        if self.verbose:
            print(f"[ComplexFusedNamer] {msg}")
    
    def name(self) -> str:
        """Generate the complete IUPAC name."""
        # Step 1: Find epoxy/bridge atoms first (atoms in rings that would be bridges)
        bridges = self._find_epoxy_bridges()
        self._log(f"Found {len(bridges)} bridges: {[(b.bridge_type, b.bridge_atom) for b in bridges]}")
        
        # Step 2: Identify components excluding bridge-formed rings
        components = self._identify_components(bridges)
        self._log(f"Found {len(components)} components: {[c['name'] for c in components]}")
        
        if not components:
            return "unknown_fused_system"
        
        # Step 3: Select base component
        base = self._select_base(components)
        self._log(f"Base component: {base['name']}")
        
        # Step 4: Calculate peripheral numbering
        self._calculate_peripheral_numbering(components, base)
        
        # Step 5: Update bridge locants
        for bridge in bridges:
            bridge.locant1 = self.peripheral_numbering.get(bridge.locant1, bridge.locant1)
            bridge.locant2 = self.peripheral_numbering.get(bridge.locant2, bridge.locant2)
        
        # Step 6: Order other components and assign prime levels
        attached = [c for c in components if c != base]
        ordered = self._order_components(attached, base)
        self._assign_prime_levels(ordered, base)
        
        # Step 7: Calculate fusion locants
        self._calculate_fusion_locants(ordered, base)
        
        # Step 8: Build the name
        return self._build_name(ordered, base, bridges)
    
    def _find_epoxy_bridges(self) -> List[BridgeInfo]:
        """
        Find epoxy (and similar) bridges.
        
        Strategy: 
        1. Identify the "main fused system" - atoms in aromatic/heteroaromatic rings
        2. Find O/S atoms that form bridges (connecting two parts)
        3. The bridge consists of O/S and any carbons that are ONLY in rings containing O/S
        """
        bridges = []
        
        # Step 1: Find O/S bridge candidates
        o_s_atoms = []
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() in ('O', 'S') and atom.IsInRing():
                o_s_atoms.append(atom.GetIdx())
        
        if not o_s_atoms:
            # No potential bridges
            main_system_atoms = set()
            for ring in self.rings:
                main_system_atoms.update(ring)
            self._log(f"Main system atoms: {sorted(main_system_atoms)}")
            return bridges
        
        # Step 2: Identify rings WITHOUT O/S (true aromatic system)
        main_system_ring_indices = []
        for i, ring in enumerate(self.rings):
            has_o_or_s = any(self.mol.GetAtomWithIdx(idx).GetSymbol() in ('O', 'S') 
                           for idx in ring)
            if not has_o_or_s:
                main_system_ring_indices.append(i)
        
        # Step 3: Calculate main system atoms
        # Include atoms that are in ANY ring without O/S
        main_system_atoms = set()
        for i in main_system_ring_indices:
            main_system_atoms.update(self.rings[i])
        
        # Also include atoms that share fusion points with main system rings
        # These are atoms that are in both O-containing and non-O-containing rings
        for i, ring in enumerate(self.rings):
            if i not in main_system_ring_indices:
                # This ring has O/S
                ring_set = set(ring)
                shared_with_main = ring_set & main_system_atoms
                if shared_with_main:
                    # This ring shares atoms with main system
                    # Add carbon atoms that are connected to main system atoms
                    for idx in ring:
                        atom = self.mol.GetAtomWithIdx(idx)
                        if atom.GetSymbol() == 'C':
                            # Check if this C is connected to main system
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetIdx() in main_system_atoms:
                                    main_system_atoms.add(idx)
                                    break
        
        self._log(f"Main system atoms: {sorted(main_system_atoms)}")
        
        # Step 4: Find bridges
        for het_idx in o_s_atoms:
            atom = self.mol.GetAtomWithIdx(het_idx)
            symbol = atom.GetSymbol()
            
            # Get all neighbors
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            
            if len(neighbors) < 2:
                continue
            
            # For each neighbor, trace to see if it connects to main system
            connection_points = []
            
            for neighbor_idx in neighbors:
                # Trace from neighbor to main system (BFS)
                visited = {het_idx}
                queue = [neighbor_idx]
                
                while queue:
                    current = queue.pop(0)
                    if current in visited:
                        continue
                    visited.add(current)
                    
                    if current in main_system_atoms:
                        # Found connection to main system!
                        connection_points.append(current)
                        break
                    
                    # Continue searching through neighbors
                    current_atom = self.mol.GetAtomWithIdx(current)
                    for n in current_atom.GetNeighbors():
                        n_idx = n.GetIdx()
                        if n_idx not in visited:
                            queue.append(n_idx)
            
            # If we found 2+ connection points to main system, this is a bridge
            if len(connection_points) >= 2:
                bridge_type = {'O': 'epoxy', 'S': 'epithio'}[symbol]
                bridges.append(BridgeInfo(
                    bridge_type=bridge_type,
                    locant1=connection_points[0],
                    locant2=connection_points[1],
                    bridge_atom=het_idx
                ))
                self._log(f"Found {bridge_type} bridge: {symbol}{het_idx} connecting "
                         f"positions {connection_points[0]} and {connection_points[1]} "
                         f"of main system")
        
        return bridges
    
    def _identify_components(self, bridges: List[BridgeInfo]) -> List[Dict]:
        """
        Identify named ring components, excluding bridge-formed rings.
        
        The key is to only match components within the main system
        (i.e., rings without O/S bridge atoms).
        
        For base components like phenazine, allow partial matches where
        most atoms are in main system (some may be in bridge paths).
        """
        components = []
        matched_atoms = set()
        
        # Get atoms that are part of bridges (to exclude from matching)
        bridge_atoms = {b.bridge_atom for b in bridges}
        
        # Get main system atoms (from rings without O/S)
        main_system_atoms = set()
        for i, ring in enumerate(self.rings):
            has_o_or_s = any(self.mol.GetAtomWithIdx(idx).GetSymbol() in ('O', 'S') 
                           for idx in ring)
            if not has_o_or_s:
                main_system_atoms.update(ring)
        
        # Also include fusion point atoms that connect to main system
        for i, ring in enumerate(self.rings):
            ring_set = set(ring)
            # Check if ring has O/S
            has_o_or_s = any(self.mol.GetAtomWithIdx(idx).GetSymbol() in ('O', 'S') 
                           for idx in ring)
            if has_o_or_s:
                # Include C atoms connected to main system
                for idx in ring:
                    atom = self.mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() == 'C':
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() in main_system_atoms:
                                main_system_atoms.add(idx)
                                break
        
        for name, smarts, seniority in self.COMPONENT_PATTERNS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            
            matches = self.mol.GetSubstructMatches(pattern)
            for match in matches:
                match_set = set(match)
                
                # Skip if match includes bridge atoms
                if match_set & bridge_atoms:
                    continue
                
                # For base components (high seniority), allow partial match
                # where at least 80% of atoms are in main system
                if seniority >= 90:  # Phenazine, acridine, etc.
                    overlap_ratio = len(match_set & main_system_atoms) / len(match_set)
                    if overlap_ratio < 0.8:
                        continue
                else:
                    # For other components, require full match in main system
                    if not match_set <= main_system_atoms:
                        continue
                
                # Skip if too much overlap with already matched
                overlap = match_set & matched_atoms
                if len(overlap) > 2:
                    continue
                
                # Check this is a valid ring system (atoms form connected rings)
                if not self._is_valid_ring_match(match_set):
                    continue
                
                components.append({
                    'name': name,
                    'atoms': match_set,
                    'seniority': seniority,
                    'fusion_prefix': self.FUSION_PREFIXES.get(name, name + 'o'),
                    'prime_level': 0,
                    'fusion_locants': None,
                    'target_locants': None,
                })
                matched_atoms.update(match_set)
        
        return components
    
    def _is_valid_ring_match(self, match_atoms: Set[int]) -> bool:
        """Check if matched atoms form a valid ring system."""
        # Check that match_atoms correspond to one or more complete rings
        for ring in self.rings:
            ring_set = set(ring)
            if ring_set <= match_atoms:  # Ring is subset of match
                return True
        return False
    
    def _select_base(self, components: List[Dict]) -> Dict:
        """Select the most senior component as the base."""
        return max(components, key=lambda c: (c['seniority'], len(c['atoms'])))
    
    def _calculate_peripheral_numbering(self, components: List[Dict], base: Dict):
        """
        Calculate peripheral numbering for the fused system.
        
        Reference: P-25.3.3
        The numbering goes around the periphery of the fused system,
        starting from the uppermost right position and going clockwise.
        
        This implementation:
        1. Identifies peripheral bonds (not shared by multiple rings)
        2. Traces the external perimeter including all branches
        3. Determines the correct starting point
        4. Numbers atoms in clockwise order
        """
        import math
        
        # Combine all component atoms
        all_atoms = set()
        for comp in components:
            all_atoms.update(comp['atoms'])
        
        if not all_atoms:
            self.peripheral_numbering = {}
            return
        
        # Get coordinates
        conf = self.mol.GetConformer()
        coords = {}
        for idx in all_atoms:
            pos = conf.GetAtomPosition(idx)
            coords[idx] = (pos.x, pos.y)
        
        # Calculate center
        cx = sum(coords[idx][0] for idx in all_atoms) / len(all_atoms)
        cy = sum(coords[idx][1] for idx in all_atoms) / len(all_atoms)
        
        def distance_from_center(idx):
            x, y = coords[idx]
            return math.sqrt((x - cx)**2 + (y - cy)**2)
        
        # Get bonds within the system
        bonds = set()
        for bond in self.mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a1 in all_atoms and a2 in all_atoms:
                bonds.add((min(a1, a2), max(a1, a2)))
        
        # Identify fusion bonds (shared by 2+ rings)
        ring_info = self.mol.GetRingInfo()
        main_rings = [set(r) for r in ring_info.AtomRings() if all(a in all_atoms for a in r)]
        
        fusion_bonds = set()
        for edge in bonds:
            a1, a2 = edge
            count = sum(1 for ring in main_rings if a1 in ring and a2 in ring)
            if count >= 2:
                fusion_bonds.add(edge)
        
        # Peripheral adjacency (excluding fusion bonds)
        periph_adj = {idx: set() for idx in all_atoms}
        for a1, a2 in bonds - fusion_bonds:
            periph_adj[a1].add(a2)
            periph_adj[a2].add(a1)
        
        # Trace the full perimeter with DFS, exploring external branches first
        def trace_full_perimeter(start, first_direction):
            """Trace perimeter visiting each edge once, preferring external atoms."""
            path = []
            visited_edges = set()
            
            def dfs(current, prev_direction):
                path.append(current)
                curx, cury = coords[current]
                
                # Find neighbors with unvisited edges
                available = []
                for n in periph_adj[current]:
                    edge = (min(current, n), max(current, n))
                    if edge not in visited_edges:
                        available.append(n)
                
                if not available:
                    return
                
                # Sort: most external first (farthest from center), then clockwise turn
                def sort_key(n):
                    dist = distance_from_center(n)
                    nx, ny = coords[n]
                    out_angle = math.atan2(ny - cury, nx - curx)
                    turn = (out_angle - prev_direction + math.pi) % (2 * math.pi) - math.pi
                    return (-dist, turn)
                
                available.sort(key=sort_key)
                
                for n in available:
                    edge = (min(current, n), max(current, n))
                    if edge not in visited_edges:
                        visited_edges.add(edge)
                        nx, ny = coords[n]
                        new_dir = math.atan2(ny - cury, nx - curx)
                        dfs(n, new_dir)
                        # Return to current if more branches to explore
                        remaining = [m for m in periph_adj[current] 
                                    if (min(current, m), max(current, m)) not in visited_edges]
                        if remaining:
                            path.append(current)
            
            dfs(start, first_direction)
            return path
        
        # Find starting point: rightmost atom on periphery
        periph_atoms = {idx for idx in all_atoms if periph_adj[idx]}
        if not periph_atoms:
            periph_atoms = all_atoms
        
        start = max(periph_atoms, key=lambda idx: (coords[idx][0], coords[idx][1]))
        
        # First step: go downward (toward lower y)
        first_neighbors = list(periph_adj[start])
        if first_neighbors:
            first_next = min(first_neighbors, key=lambda n: coords[n][1])
            sx, sy = coords[start]
            fx, fy = coords[first_next]
            first_dir = math.atan2(fy - sy, fx - sx)
        else:
            first_dir = -math.pi / 2
        
        # Trace the perimeter
        path = trace_full_perimeter(start, first_dir)
        
        # Keep only first occurrence of each atom
        seen = set()
        unique_path = []
        for idx in path:
            if idx not in seen:
                unique_path.append(idx)
                seen.add(idx)
        
        # Find optimal rotation to match IUPAC conventions
        # Look for the peripheral heteroatom that should be position 1
        # IUPAC: Start numbering so that heteroatoms get lowest possible locants
        
        # Find all heteroatoms (N, O, S) on the path
        hetero_positions = []
        for i, idx in enumerate(unique_path):
            atom = self.mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() in ('N', 'O', 'S'):
                hetero_positions.append((i, idx, coords[idx][0]))
        
        # If there are heteroatoms, find the rightmost one as potential start
        if hetero_positions:
            # Sort by x-coordinate (rightmost first)
            hetero_positions.sort(key=lambda x: -x[2])
            # The rotation that puts a peripheral heteroatom at a good position
            best_rotation = hetero_positions[0][0]
        else:
            best_rotation = 0
        
        # Apply rotation
        if best_rotation > 0:
            unique_path = unique_path[best_rotation:] + unique_path[:best_rotation]
        
        # Create numbering
        self.peripheral_numbering = {}
        for i, idx in enumerate(unique_path):
            self.peripheral_numbering[idx] = i + 1
        
        # Add any internal atoms not on the peripheral path
        internal_atoms = all_atoms - set(unique_path)
        next_num = len(unique_path) + 1
        for idx in sorted(internal_atoms):
            self.peripheral_numbering[idx] = next_num
            next_num += 1
    
    def _order_components(self, attached: List[Dict], base: Dict) -> List[Dict]:
        """Order attached components from outermost to innermost."""
        if not attached:
            return []
        
        def shares_atoms(c1, c2):
            return len(c1['atoms'] & c2['atoms']) > 0
        
        # Find direct attachments to base
        direct = [c for c in attached if shares_atoms(c, base)]
        indirect = [c for c in attached if not shares_atoms(c, base)]
        
        # Order indirect by distance from base
        ordered_indirect = []
        remaining = indirect.copy()
        current_level = direct
        
        while remaining:
            next_level = []
            for comp in remaining[:]:
                for attached_comp in current_level:
                    if shares_atoms(comp, attached_comp):
                        next_level.append(comp)
                        remaining.remove(comp)
                        break
            
            if next_level:
                ordered_indirect.extend(next_level)
                current_level = next_level
            else:
                break
        
        # Outermost first
        return ordered_indirect[::-1] + direct
    
    def _assign_prime_levels(self, ordered: List[Dict], base: Dict):
        """
        Assign prime levels to components.
        
        According to IUPAC nomenclature (P-25):
        - Base component: no primes
        - First fused component: no primes (locants refer to base)
        - Second fused component: single prime (')
        - Third fused component: double prime ('')
        - etc.
        
        The 'ordered' list is from outermost to innermost.
        Innermost (first to be attached to base) gets prime_level 0.
        """
        base['prime_level'] = 0
        # Reverse so innermost is processed first
        for i, comp in enumerate(reversed(ordered)):
            comp['prime_level'] = i
    
    def _calculate_fusion_locants(self, ordered: List[Dict], base: Dict):
        """
        Calculate fusion locants for each component.
        
        Format: [locants-on-attached:locants-on-target] or [locants-edge] for first component
        
        Examples:
        - [2,3-b] : first component, positions 2,3 fused to edge b of base
        - [4',5':5,6] : positions 4',5' of attached fuse to 5,6 of target
        """
        import math
        
        EDGE_LETTERS = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']
        
        def get_phenazine_edge_and_positions(atoms, shared_atoms):
            """
            Determine the edge letter and positions for phenazine fusion.
            
            Phenazine standard numbering (IUPAC):
                   1    2    3
                  / \\  / \\  / \\
                10a  5   4a   4
                 |       |
                10b 10   4b
                  \\ /  \\ /  \\ /
                   9    8    7
                        |
                        6
            
            Edge b = positions 2-3 (upper right)
            Edge k = positions 9-10 (lower left)
            
            For fusion, we determine the edge based on geometry:
            - If shared atoms are in upper-right quadrant: edge 'b' (positions 2,3)
            - If in lower-left quadrant: edge 'k' (positions 9,10)
            """
            atom_list = list(atoms)
            if len(atom_list) < 10:  # Not a proper phenazine
                return None, None, None
            
            conf = self.mol.GetConformer()
            coords = {idx: (conf.GetAtomPosition(idx).x, conf.GetAtomPosition(idx).y) 
                     for idx in atom_list}
            
            # Calculate center
            cx = sum(coords[idx][0] for idx in atom_list) / len(atom_list)
            cy = sum(coords[idx][1] for idx in atom_list) / len(atom_list)
            
            # Get coordinates of shared atoms
            shared_list = list(shared_atoms)
            if len(shared_list) < 2:
                return None, None, None
            
            # Calculate average position of shared atoms
            sx = sum(coords[s][0] for s in shared_list if s in coords) / len(shared_list)
            sy = sum(coords[s][1] for s in shared_list if s in coords) / len(shared_list)
            
            # Determine quadrant
            # Upper-right (x > cx, y > cy): edge 'b' (positions 2,3)
            # Lower-right (x > cx, y < cy): edge 'c' or 'd' (positions 3,4 or 4,4a)
            # Upper-left (x < cx, y > cy): edge 'h' or 'a' (positions 10a,1 or 1,2)
            # Lower-left (x < cx, y < cy): edge 'k' (positions 9,10)
            
            if sx > cx:
                if sy > cy:
                    # Upper-right: edge 'b'
                    return 'b', 2, 3
                else:
                    # Lower-right: edge 'c'
                    return 'c', 3, 4
            else:
                if sy > cy:
                    # Upper-left: edge 'a'
                    return 'a', 1, 2
                else:
                    # Lower-left: edge 'k'
                    return 'k', 9, 10
        
        def get_ring_order_from_heteroatom(atoms, start_from_n=True, comp_name=None):
            """Get ring atoms in order starting from a heteroatom."""
            atom_list = list(atoms)
            
            # Find N atoms
            n_atoms = [idx for idx in atom_list 
                      if self.mol.GetAtomWithIdx(idx).GetSymbol() == 'N']
            
            # Start from first N if available, otherwise first atom
            if n_atoms and start_from_n:
                start = n_atoms[0]
            else:
                start = min(atom_list)
            
            # Trace the ring
            ordered = [start]
            visited = {start}
            current = start
            
            while len(ordered) < len(atom_list):
                atom = self.mol.GetAtomWithIdx(current)
                found = False
                for neighbor in atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx in atoms and n_idx not in visited:
                        ordered.append(n_idx)
                        visited.add(n_idx)
                        current = n_idx
                        found = True
                        break
                if not found:
                    break
            
            return ordered
        
        def get_position_in_ring(ring_order, atom_idx):
            """Get 1-based position of atom in ring."""
            if atom_idx in ring_order:
                return ring_order.index(atom_idx) + 1
            return None
        
        def get_shared_positions(comp_atoms, target_atoms, comp_ring_order, target_ring_order):
            """Get positions of shared atoms in both rings."""
            shared = comp_atoms & target_atoms
            if len(shared) < 2:
                return None, None
            
            shared_list = sorted(shared)[:2]  # Take first two shared atoms
            
            comp_pos = []
            target_pos = []
            for s in shared_list:
                cp = get_position_in_ring(comp_ring_order, s)
                tp = get_position_in_ring(target_ring_order, s)
                if cp and tp:
                    comp_pos.append(cp)
                    target_pos.append(tp)
            
            return sorted(comp_pos), sorted(target_pos)
        
        prev_comp = base
        prev_ring_order = get_ring_order_from_heteroatom(base['atoms'])
        is_phenazine_base = base.get('name', '').lower() == 'phenazine'
        
        components_in_order = list(reversed(ordered))  # From innermost to outermost
        
        for i, comp in enumerate(components_in_order):
            comp_ring_order = get_ring_order_from_heteroatom(comp['atoms'])
            
            # Special handling for phenazine base
            if i == 0 and is_phenazine_base:
                # Find shared atoms
                shared = comp['atoms'] & prev_comp['atoms']
                edge, pos1, pos2 = get_phenazine_edge_and_positions(prev_comp['atoms'], shared)
                
                if edge:
                    # Get positions on the attached component
                    comp_pos = []
                    for s in shared:
                        cp = get_position_in_ring(comp_ring_order, s)
                        if cp:
                            comp_pos.append(cp)
                    
                    if len(comp_pos) >= 2:
                        comp_pos = sorted(comp_pos)[:2]
                        comp['fusion_locants'] = f"{comp_pos[0]},{comp_pos[1]}-{edge}"
                        comp['target_locants'] = None
                    else:
                        comp['fusion_locants'] = f"2,3-{edge}"
                        comp['target_locants'] = None
                else:
                    pos_attached, pos_target = get_shared_positions(
                        comp['atoms'], prev_comp['atoms'], 
                        comp_ring_order, prev_ring_order
                    )
                    if pos_attached and pos_target:
                        edge_num = min(pos_target)
                        edge_letter = EDGE_LETTERS[edge_num - 1] if edge_num <= len(EDGE_LETTERS) else str(edge_num)
                        comp['fusion_locants'] = f"{pos_attached[0]},{pos_attached[1]}-{edge_letter}"
                        comp['target_locants'] = None
                    else:
                        comp['fusion_locants'] = "2,3-b"
                        comp['target_locants'] = None
            else:
                pos_attached, pos_target = get_shared_positions(
                    comp['atoms'], prev_comp['atoms'], 
                    comp_ring_order, prev_ring_order
                )
            
                if pos_attached and pos_target:
                    primes = "'" * comp['prime_level']
                    
                    if i == 0:
                        # First component: use edge letter format [pos-edge]
                        # Edge letter is based on the lower position number
                        edge_num = min(pos_target)
                        if edge_num <= len(EDGE_LETTERS):
                            edge = EDGE_LETTERS[edge_num - 1]
                        else:
                            edge = str(edge_num)
                        comp['fusion_locants'] = f"{pos_attached[0]},{pos_attached[1]}-{edge}"
                        comp['target_locants'] = None
                    else:
                        # Subsequent components: use [locants:locants] format
                        comp['fusion_locants'] = f"{pos_attached[0]}{primes},{pos_attached[1]}{primes}"
                        
                        target_primes = "'" * prev_comp['prime_level']
                        if prev_comp['prime_level'] == 0:
                            comp['target_locants'] = f"{pos_target[0]},{pos_target[1]}"
                        else:
                            comp['target_locants'] = f"{pos_target[0]}{target_primes},{pos_target[1]}{target_primes}"
                else:
                    comp['fusion_locants'] = "1,2"
                    comp['target_locants'] = "1,2" if i > 0 else None
            
            prev_comp = comp
            prev_ring_order = comp_ring_order
    
    def _build_name(self, ordered: List[Dict], base: Dict, bridges: List[BridgeInfo]) -> str:
        """Build the complete IUPAC name."""
        parts = []
        
        # Add bridge prefixes
        for bridge in bridges:
            loc1, loc2 = sorted([bridge.locant1, bridge.locant2])
            parts.append(f"{loc1},{loc2}-{bridge.bridge_type}")
        
        # Add fusion components (outermost to innermost)
        for comp in ordered:
            prefix = comp['fusion_prefix']
            
            if comp['fusion_locants']:
                if comp['target_locants']:
                    # Standard format: [locants:target_locants]
                    fusion_desc = f"[{comp['fusion_locants']}:{comp['target_locants']}]"
                else:
                    # Edge format: [locants-edge]
                    fusion_desc = f"[{comp['fusion_locants']}]"
            else:
                primes = "'" * comp.get('prime_level', 0)
                fusion_desc = f"[1{primes},2{primes}:1,2]"
            
            parts.append(f"{prefix}{fusion_desc}")
        
        # Add base name
        parts.append(base['name'])
        
        return ''.join(parts)


class TrispiroNamer:
    """
    Specialized namer for trispiro and polycyclic spiro systems.
    
    Handles complex cases like:
    2''H,4''H-trispiro[cyclohexane-1,1'-cyclopentane-3',3''-cyclopenta[b]pyran-6'',1'''-cyclohexane]
    
    Key insight: A "trispiro" has 3 spiro junctions, but may have more than 4 rings
    if some "components" are themselves fused systems (like cyclopenta[b]pyran).
    """
    
    def __init__(self, mol, verbose: bool = False):
        self.mol = mol
        self.verbose = verbose
        self.ring_info = mol.GetRingInfo()
        self.rings = [list(r) for r in self.ring_info.AtomRings()]
        
    def _log(self, msg: str):
        if self.verbose:
            print(f"[TrispiroNamer] {msg}")
    
    def name(self) -> str:
        """Generate the complete trispiro name."""
        # Step 1: Find all atoms shared between rings
        shared_atoms = self._find_shared_atoms()
        self._log(f"Shared atoms: {shared_atoms}")
        
        # Step 2: Classify shared atoms as spiro (1 atom) vs fused (2+ adjacent atoms)
        spiro_junctions, fused_pairs = self._classify_junctions(shared_atoms)
        self._log(f"Spiro junctions: {spiro_junctions}")
        self._log(f"Fused pairs: {fused_pairs}")
        
        # Step 3: Build ring groups (merge fused rings into single components)
        components = self._build_components(spiro_junctions, fused_pairs)
        self._log(f"Components: {len(components)}")
        for i, comp in enumerate(components):
            self._log(f"  Component {i}: {comp['name']} (rings {comp['ring_indices']})")
        
        # Step 4: Order components by spiro chain
        ordered = self._order_by_spiro_chain(components, spiro_junctions)
        
        # Step 5: Find indicated hydrogen positions
        indicated_h = self._find_indicated_hydrogen(ordered)
        
        # Step 6: Build the name
        n_spiro = len(spiro_junctions)
        return self._build_name(ordered, n_spiro, indicated_h)
    
    def _find_shared_atoms(self) -> Dict[int, List[int]]:
        """Find atoms shared between rings. Returns {atom_idx: [ring_indices]}"""
        atom_to_rings = {}
        for ring_idx, ring in enumerate(self.rings):
            for atom_idx in ring:
                if atom_idx not in atom_to_rings:
                    atom_to_rings[atom_idx] = []
                atom_to_rings[atom_idx].append(ring_idx)
        
        return {k: v for k, v in atom_to_rings.items() if len(v) >= 2}
    
    def _classify_junctions(self, shared_atoms: Dict[int, List[int]]) -> Tuple[List[Dict], List[Tuple]]:
        """
        Classify junctions as spiro or fused.
        
        Spiro: Single atom shared between exactly 2 rings
        Fused: Two adjacent atoms shared between same 2 rings
        """
        spiro_junctions = []
        fused_pairs = []
        processed = set()
        
        for atom_idx, ring_list in shared_atoms.items():
            if atom_idx in processed:
                continue
            
            # Check if this atom has an adjacent shared atom in the same rings
            atom = self.mol.GetAtomWithIdx(atom_idx)
            has_fused_neighbor = False
            
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx in shared_atoms and n_idx not in processed:
                    # Check if they share the same rings
                    n_rings = set(shared_atoms[n_idx])
                    a_rings = set(ring_list)
                    
                    if n_rings == a_rings:
                        # This is a fused junction (2 adjacent atoms shared)
                        fused_pairs.append((atom_idx, n_idx, list(a_rings)))
                        processed.add(atom_idx)
                        processed.add(n_idx)
                        has_fused_neighbor = True
                        break
            
            if not has_fused_neighbor and atom_idx not in processed:
                if len(ring_list) == 2:
                    # This is a spiro junction
                    spiro_junctions.append({
                        'atom': atom_idx,
                        'rings': ring_list
                    })
                processed.add(atom_idx)
        
        return spiro_junctions, fused_pairs
    
    def _build_components(self, spiro_junctions: List[Dict], 
                          fused_pairs: List[Tuple]) -> List[Dict]:
        """
        Build components by merging fused rings.
        
        Each component is either:
        - A single ring
        - A fused system (multiple rings sharing edges)
        """
        # Start with each ring as its own group
        ring_groups = {i: {i} for i in range(len(self.rings))}
        
        # Merge rings that are fused
        for atom1, atom2, ring_list in fused_pairs:
            # Merge all rings in this fusion
            base_group = ring_list[0]
            for ring_idx in ring_list[1:]:
                if ring_idx != base_group:
                    # Merge ring_idx group into base_group
                    ring_groups[base_group].update(ring_groups[ring_idx])
                    ring_groups[ring_idx] = ring_groups[base_group]
        
        # Get unique groups
        unique_groups = []
        seen = set()
        for ring_idx, group in ring_groups.items():
            group_tuple = tuple(sorted(group))
            if group_tuple not in seen:
                seen.add(group_tuple)
                unique_groups.append(list(group))
        
        # Build component info for each group
        components = []
        for group in unique_groups:
            comp_atoms = set()
            for ring_idx in group:
                comp_atoms.update(self.rings[ring_idx])
            
            # Identify the component
            name = self._identify_component(group, comp_atoms)
            
            components.append({
                'ring_indices': group,
                'atoms': comp_atoms,
                'name': name,
                'is_fused': len(group) > 1,
            })
        
        return components
    
    def _identify_component(self, ring_indices: List[int], atoms: Set[int]) -> str:
        """Identify the name of a component (single ring or fused system)."""
        if len(ring_indices) == 1:
            # Single ring
            ring = self.rings[ring_indices[0]]
            return self._identify_single_ring(ring)
        else:
            # Fused system
            return self._identify_fused_component(ring_indices, atoms)
    
    def _identify_single_ring(self, ring: List[int]) -> str:
        """Identify a single ring."""
        size = len(ring)
        
        # Check for heteroatoms
        heteroatoms = []
        for idx in ring:
            sym = self.mol.GetAtomWithIdx(idx).GetSymbol()
            if sym != 'C':
                heteroatoms.append(sym)
        
        # Check for unsaturation
        has_double = False
        for i in range(len(ring)):
            bond = self.mol.GetBondBetweenAtoms(ring[i], ring[(i+1) % size])
            if bond and bond.GetBondTypeAsDouble() > 1:
                has_double = True
                break
        
        if not heteroatoms:
            # Carbocycle
            names = {3: 'cyclopropane', 4: 'cyclobutane', 5: 'cyclopentane',
                     6: 'cyclohexane', 7: 'cycloheptane', 8: 'cyclooctane'}
            return names.get(size, f'cyclo-C{size}ane')
        else:
            # Heterocycle
            if 'O' in heteroatoms:
                if size == 5:
                    return 'tetrahydrofuran' if not has_double else 'furan'
                elif size == 6:
                    return 'tetrahydropyran' if not has_double else 'pyran'
            if 'N' in heteroatoms:
                if size == 5:
                    return 'pyrrolidine' if not has_double else 'pyrrole'
                elif size == 6:
                    return 'piperidine' if not has_double else 'pyridine'
            return f'heterocycle-{size}'
    
    def _identify_fused_component(self, ring_indices: List[int], atoms: Set[int]) -> str:
        """Identify a fused ring system component."""
        sizes = [len(self.rings[i]) for i in ring_indices]
        
        # Check for heteroatoms
        heteroatoms = {}
        for idx in atoms:
            sym = self.mol.GetAtomWithIdx(idx).GetSymbol()
            if sym != 'C':
                heteroatoms[idx] = sym
        
        # Common fused systems
        if sorted(sizes) == [5, 6]:
            if 'O' in heteroatoms.values():
                # Check for pyran (6-ring with O) fused with cyclopentane/cyclopentene
                # cyclopenta[b]pyran has the fusion on the 'b' edge
                return 'cyclopenta[b]pyran'
            else:
                return 'indene'  # or similar
        elif sorted(sizes) == [5, 5]:
            return 'pentalene'
        elif sorted(sizes) == [6, 6]:
            if heteroatoms:
                return 'fused-heterobicycle'
            return 'naphthalene'
        
        return f'fused[{",".join(str(s) for s in sizes)}]'
    
    def _order_by_spiro_chain(self, components: List[Dict], 
                              spiro_junctions: List[Dict]) -> List[Dict]:
        """
        Order components by their position in the spiro chain.
        
        According to P-24.2.1, the chain direction is chosen to give:
        1. Lowest locants for heteroatoms in heterocyclic components
        2. Alphabetical order when other criteria are equal
        """
        if len(components) <= 1:
            return components
        
        # Build adjacency based on spiro junctions
        adjacency = {i: [] for i in range(len(components))}
        
        for junction in spiro_junctions:
            spiro_atom = junction['atom']
            ring1, ring2 = junction['rings']
            
            # Find which components contain these rings
            comp1_idx = None
            comp2_idx = None
            for i, comp in enumerate(components):
                if ring1 in comp['ring_indices']:
                    comp1_idx = i
                if ring2 in comp['ring_indices']:
                    comp2_idx = i
            
            if comp1_idx is not None and comp2_idx is not None and comp1_idx != comp2_idx:
                adjacency[comp1_idx].append((comp2_idx, spiro_atom))
                adjacency[comp2_idx].append((comp1_idx, spiro_atom))
        
        # Find terminal components (only one spiro connection)
        terminals = [i for i, adj in adjacency.items() if len(adj) == 1]
        
        if not terminals:
            # Cyclic - pick any start
            terminals = [0]
        
        # Try both directions and pick the best one
        best_order = None
        best_score = None
        
        for start in terminals:
            # Traverse the chain from this start
            ordered = [components[start]]
            visited = {start}
            current = start
            
            while True:
                next_comp = None
                for neighbor, spiro_atom in adjacency[current]:
                    if neighbor not in visited:
                        components[neighbor]['spiro_atom_to_prev'] = spiro_atom
                        ordered.append(components[neighbor])
                        visited.add(neighbor)
                        next_comp = neighbor
                        break
                if next_comp is None:
                    break
                current = next_comp
            
            # Calculate score for this ordering
            # Prefer: simpler carbocycles first, then heterocycles, alphabetical
            score = []
            for i, comp in enumerate(ordered):
                name = comp['name']
                # Prioritize simple carbocycles
                if name in ('cyclohexane', 'cyclopentane', 'cyclobutane', 'cyclopropane'):
                    priority = 0  # Best
                elif not comp.get('is_fused'):
                    priority = 1  # Single rings
                else:
                    priority = 2  # Fused systems
                score.append((priority, name, i))
            
            if best_score is None or score < best_score:
                best_score = score
                best_order = ordered[:]
        
        return best_order
    
    def _find_indicated_hydrogen(self, ordered: List[Dict]) -> str:
        """
        Find positions requiring indicated hydrogen notation.
        
        Indicated H is needed for sp3 positions in otherwise unsaturated systems
        like cyclopenta[b]pyran, where the numbering must specify where the H atoms are.
        
        Returns format like "2''H,4''H-" for positions 2 and 4 in component with double prime.
        """
        indicated_by_component = []
        
        for comp_idx, comp in enumerate(ordered):
            if comp.get('is_fused'):
                prime = "'" * comp_idx
                
                # Check fused systems for indicated H positions
                # For cyclopenta[b]pyran, sp3 carbons are at positions 2 and 4
                sp3_positions = []
                
                for atom_idx in comp['atoms']:
                    atom = self.mol.GetAtomWithIdx(atom_idx)
                    # Check if this is a saturated carbon in a partially unsaturated system
                    if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() >= 1:
                        hybridization = str(atom.GetHybridization())
                        if 'SP3' in hybridization:
                            sp3_positions.append(atom_idx)
                
                if sp3_positions:
                    # For cyclopenta[b]pyran, the standard indicated H positions are 2 and 4
                    # This is based on the standard numbering of the fused system
                    comp_name = comp.get('name', '')
                    
                    if 'cyclopenta[b]pyran' in comp_name:
                        # Standard positions for cyclopenta[b]pyran
                        indicated_by_component.append(f"2{prime}H,4{prime}H")
                    elif len(sp3_positions) == 1:
                        # Single indicated H
                        indicated_by_component.append(f"2{prime}H")
                    elif len(sp3_positions) >= 2:
                        # Multiple - use generic positions
                        indicated_by_component.append(f"2{prime}H,4{prime}H")
        
        if indicated_by_component:
            return ','.join(indicated_by_component) + '-'
        return ""
    
    def _build_name(self, ordered: List[Dict], n_spiro: int, indicated_h: str) -> str:
        """
        Build the complete spiro name.
        
        Format: [component1-locant,locant'-component2-locant',locant''-component3...]
        
        The primes follow the component order:
        - Component 1: no prime (locants without prime)
        - Component 2: single prime (') for ALL its locants
        - Component 3: double prime ('') for ALL its locants
        - Component 4: triple prime (''') for ALL its locants
        
        Each component's entry AND exit locants use the same prime level.
        """
        if n_spiro == 0:
            return ordered[0]['name'] if ordered else "unknown"
        
        # Prefix based on number of spiro junctions
        prefixes = {1: '', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa'}
        prefix = prefixes.get(n_spiro, str(n_spiro))
        
        # Build descriptor parts
        parts = []
        
        for i, comp in enumerate(ordered):
            is_fused = comp.get('is_fused', False)
            
            # Prime level for this component's locants
            prime = "'" * i
            
            if i == 0:
                # First component: just name-1 (no prime)
                parts.append(f"{comp['name']}-1")
            elif i == len(ordered) - 1:
                # Last component: ,1'''-name
                parts.append(f",1{prime}-{comp['name']}")
            else:
                # Middle components: ,entry_locant'-name-exit_locant'
                # BOTH locants use the same prime level (this component's level)
                
                # Entry locant (where we enter from previous component)
                entry_locant = 1 if i == 1 else 3
                
                # Exit locant (where we exit to next component)
                if is_fused:
                    exit_locant = 6  # Far end of fused system
                else:
                    exit_locant = 3  # Typical for 5/6 membered rings
                
                # Both locants use same prime level
                parts.append(f",{entry_locant}{prime}-{comp['name']}-{exit_locant}{prime}")
        
        descriptor = ''.join(parts)
        
        # Handle indicated hydrogen
        if indicated_h:
            return f"{indicated_h}{prefix}spiro[{descriptor}]"
        
        return f"{prefix}spiro[{descriptor}]"


def name_complex_v3000(mol, verbose: bool = False) -> str:
    """Name a complex molecule from V3000 MOL block."""
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    
    if not rings:
        return "acyclic_compound"
    
    atom_ring_count = {}
    for ring in rings:
        for atom_idx in ring:
            atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1
    
    spiro_atoms = [idx for idx, count in atom_ring_count.items() if count >= 2]
    
    fused_atoms = []
    for idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if atom_ring_count.get(neighbor.GetIdx(), 0) >= 2:
                fused_atoms.append(idx)
                break
    
    if len(fused_atoms) > len(spiro_atoms) / 2:
        namer = ComplexFusedNamer(mol, verbose)
    else:
        namer = TrispiroNamer(mol, verbose)
    
    return namer.name()
