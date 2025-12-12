"""
P-24: Spiro Ring Systems
========================

Comprehensive rules for naming spiro compounds.

Reference: IUPAC Blue Book 2013, P-24

Spiro compounds have two or more rings sharing exactly one atom (the spiro atom).
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass


@dataclass
class SpiroSystem:
    """Represents a spiro ring system."""
    spiro_atoms: List[int]          # Atoms shared between rings
    rings: List[List[int]]          # Ring atom lists
    ring_sizes: List[int]           # Size of each ring
    heteroatoms: Dict[int, str]     # Position -> heteroatom symbol
    is_monospiro: bool              # Single spiro atom
    is_dispiro: bool                # Two spiro atoms
    is_trispiro: bool               # Three spiro atoms
    

# =============================================================================
# P-24.1: MONOSPIRO COMPOUNDS
# =============================================================================

def name_monospiro(ring_sizes: List[int], parent_name: str = None) -> str:
    """
    Name a monospiro compound (single spiro atom).
    
    Format: spiro[x.y]parent
    
    Where x and y are the number of atoms in each ring, not counting
    the spiro atom.
    
    Example: spiro[4.5]decane (5-membered + 6-membered sharing 1 atom)
    
    Args:
        ring_sizes: List of ring sizes
        parent_name: Optional parent hydride name
        
    Returns:
        Spiro name
    """
    if len(ring_sizes) != 2:
        return None
    
    # Ring descriptors are ring_size - 1 (excluding spiro atom)
    descriptors = sorted([s - 1 for s in ring_sizes])
    descriptor_str = '.'.join(str(d) for d in descriptors)
    
    # Total atoms = sum of ring sizes - 1 (spiro atom counted once)
    total_atoms = sum(ring_sizes) - 1
    
    if parent_name is None:
        parent_name = get_alkane_name(total_atoms)
    
    return f"spiro[{descriptor_str}]{parent_name}"


def get_alkane_name(n_carbons: int) -> str:
    """Get alkane name for n carbons."""
    names = {
        1: 'methane', 2: 'ethane', 3: 'propane', 4: 'butane',
        5: 'pentane', 6: 'hexane', 7: 'heptane', 8: 'octane',
        9: 'nonane', 10: 'decane', 11: 'undecane', 12: 'dodecane',
        13: 'tridecane', 14: 'tetradecane', 15: 'pentadecane',
        16: 'hexadecane', 17: 'heptadecane', 18: 'octadecane',
        19: 'nonadecane', 20: 'icosane', 21: 'henicosane',
        22: 'docosane', 23: 'tricosane'
    }
    return names.get(n_carbons, f"C{n_carbons}ane")


# =============================================================================
# P-24.2: POLYSPIRO COMPOUNDS (DISPIRO, TRISPIRO, etc.)
# =============================================================================

def name_polyspiro(spiro_atoms: List[int], rings: List[List[int]], 
                   mol=None) -> str:
    """
    Name a polyspiro compound (multiple spiro atoms).
    
    Format: dispiro[...]name or trispiro[...]name
    
    The descriptor shows the path through all rings.
    
    Example: dispiro[2.0.2.1]heptane
    Example: trispiro[2.0.2.1.2.1]undecane
    
    Args:
        spiro_atoms: List of spiro atom indices
        rings: List of ring atom lists
        mol: Optional RDKit molecule for heteroatom detection
        
    Returns:
        Polyspiro name
    """
    n_spiro = len(spiro_atoms)
    
    if n_spiro == 1:
        ring_sizes = [len(r) for r in rings]
        return name_monospiro(ring_sizes)
    
    # Determine prefix
    prefixes = {2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa'}
    prefix = prefixes.get(n_spiro, f'{n_spiro}')
    
    # Build descriptor
    # The descriptor is a sequence of numbers showing atoms between spiro centers
    descriptor = build_polyspiro_descriptor(spiro_atoms, rings)
    
    # Calculate total atoms
    all_atoms = set()
    for ring in rings:
        all_atoms.update(ring)
    total_atoms = len(all_atoms)
    
    parent = get_alkane_name(total_atoms)
    
    return f"{prefix}spiro[{descriptor}]{parent}"


def build_polyspiro_descriptor(spiro_atoms: List[int], 
                               rings: List[List[int]]) -> str:
    """
    Build the numeric descriptor for a polyspiro compound.
    
    The descriptor format follows a path through all rings, showing
    the number of atoms between adjacent spiro atoms.
    """
    # Simplified: return ring sizes minus overlaps
    # Full implementation would trace the path properly
    
    parts = []
    spiro_set = set(spiro_atoms)
    
    for ring in rings:
        # Count non-spiro atoms in this ring
        non_spiro = len([a for a in ring if a not in spiro_set])
        parts.append(str(non_spiro))
    
    return '.'.join(parts)


# =============================================================================
# P-24.3: SPIRO HETEROCYCLES  
# =============================================================================

def name_spiro_heterocycle(rings: List[List[int]], mol,
                          heteroatom_positions: Dict[int, str]) -> str:
    """
    Name a spiro compound containing heteroatoms.
    
    Uses von Baeyer or Hantzsch-Widman names for the component rings.
    
    Example: 2-oxaspiro[4.5]decane
    
    Args:
        rings: List of ring atom lists
        mol: RDKit molecule
        heteroatom_positions: Dict mapping position to heteroatom symbol
        
    Returns:
        Heterocyclic spiro name
    """
    ring_sizes = [len(r) for r in rings]
    
    # Check for named heterocyclic components
    components = identify_spiro_components(rings, mol)
    
    if components:
        return name_spiro_with_components(components, mol)
    
    # Fallback to oxa/aza prefixes
    base_name = name_monospiro(ring_sizes) if len(ring_sizes) == 2 else \
                name_polyspiro([], rings, mol)
    
    # Add heteroatom prefixes
    prefixes = build_heteroatom_prefixes(heteroatom_positions)
    
    return f"{prefixes}{base_name}"


def identify_spiro_components(rings: List[List[int]], mol) -> List[Dict]:
    """
    Identify named ring components in a spiro system.
    
    Returns list of component info dicts.
    """
    from analyzers.component_matcher import ComponentMatcher
    
    components = []
    matcher = ComponentMatcher()
    
    for ring in rings:
        ring_set = set(ring)
        match = matcher.identify_ring(ring, mol)
        if match:
            components.append({
                'name': match,
                'atoms': ring_set,
                'size': len(ring)
            })
    
    return components


def name_spiro_with_components(components: List[Dict], mol) -> str:
    """
    Build spiro name using identified ring components.
    
    Format: spiro[componentA-locantA,locantB'-componentB]
    
    Example: spiro[cyclohexane-1,1'-cyclopentane]
    """
    if len(components) == 2:
        # Simple dispiro with named components
        c1, c2 = components[0], components[1]
        
        # Determine spiro positions (simplified)
        # The locant indicates which position is the spiro junction
        return f"spiro[{c1['name']}-1,1'-{c2['name']}]"
    
    elif len(components) > 2:
        # Complex polyspiro
        return build_complex_spiro_name(components, mol)
    
    return "spiro-unknown"


def build_complex_spiro_name(components: List[Dict], mol) -> str:
    """
    Build name for complex trispiro or higher systems.
    
    Format for trispiro:
    trispiro[componentA-1,1'-componentB-3',3''-componentC-1'',1'''-componentD]
    
    The primes indicate the numbering within each component ring.
    """
    n_components = len(components)
    
    # Get prefix
    prefixes = {2: '', 3: 'tri', 4: 'tetra', 5: 'penta'}
    prefix = prefixes.get(n_components - 1, f'{n_components-1}')
    
    parts = []
    prime_count = 0
    
    for i, comp in enumerate(components):
        name = comp['name']
        
        if i == 0:
            # First component
            parts.append(f"{name}-1")
        elif i == len(components) - 1:
            # Last component
            primes = "'" * prime_count
            parts.append(f",1{primes}-{name}")
        else:
            # Middle components
            in_primes = "'" * (prime_count)
            prime_count += 1
            out_primes = "'" * prime_count
            # Determine the locant within this ring (simplified)
            inner_locant = 3  # Would need proper analysis
            parts.append(f",{inner_locant}{in_primes}-{name}-{inner_locant}{out_primes}")
    
    return f"{prefix}spiro[{''.join(parts)}]"


def build_heteroatom_prefixes(positions: Dict[int, str]) -> str:
    """
    Build heteroatom replacement prefixes (oxa, aza, thia).
    
    Reference: P-24.3.1
    """
    prefix_names = {
        'O': 'oxa',
        'N': 'aza', 
        'S': 'thia',
        'Se': 'selena',
        'Te': 'tellura',
        'P': 'phospha',
        'As': 'arsa',
        'Si': 'sila',
        'Ge': 'germa',
        'B': 'bora',
    }
    
    # Group by heteroatom type
    by_type = {}
    for pos, symbol in positions.items():
        if symbol not in by_type:
            by_type[symbol] = []
        by_type[symbol].append(pos)
    
    # Build prefixes in order O > S > Se > Te > N > P > As > Si > Ge > B
    order = ['O', 'S', 'Se', 'Te', 'N', 'P', 'As', 'Si', 'Ge', 'B']
    
    parts = []
    for symbol in order:
        if symbol in by_type:
            locs = sorted(by_type[symbol])
            loc_str = ','.join(str(l) for l in locs)
            prefix = prefix_names.get(symbol, symbol.lower() + 'a')
            
            if len(locs) > 1:
                mult = {2: 'di', 3: 'tri', 4: 'tetra'}.get(len(locs), '')
                parts.append(f"{loc_str}-{mult}{prefix}")
            else:
                parts.append(f"{loc_str}-{prefix}")
    
    return ''.join(parts)


# =============================================================================
# P-24.4: SPIRO NUMBERING
# =============================================================================

def number_spiro_system(spiro_atoms: List[int], rings: List[List[int]],
                       mol=None) -> Dict[int, int]:
    """
    Number atoms in a spiro system.
    
    Rules (P-24.4):
    1. Start at terminal ring (ring bonded to only one other)
    2. Number continuously through all rings
    3. Spiro atoms get lowest possible numbers
    4. Heteroatoms get lowest possible numbers within constraints
    
    Args:
        spiro_atoms: List of spiro atom indices
        rings: List of ring atom lists
        mol: Optional RDKit molecule
        
    Returns:
        Dict mapping atom index to position number
    """
    numbering = {}
    current_num = 1
    visited = set()
    
    # Find terminal rings (connected to only one other ring via spiro)
    ring_connections = analyze_ring_connections(rings, spiro_atoms)
    
    # Start from terminal ring
    terminal_rings = [i for i, conns in enumerate(ring_connections) if len(conns) == 1]
    
    if not terminal_rings:
        # All rings connected in a cycle - pick smallest
        start_ring_idx = 0
    else:
        start_ring_idx = terminal_rings[0]
    
    # Number starting from first atom after spiro in terminal ring
    # ... (simplified implementation)
    
    for ring in rings:
        for atom in ring:
            if atom not in visited:
                numbering[atom] = current_num
                current_num += 1
                visited.add(atom)
    
    return numbering


def analyze_ring_connections(rings: List[List[int]], 
                            spiro_atoms: List[int]) -> List[Set[int]]:
    """
    Determine which rings are connected via spiro atoms.
    
    Returns list of sets, where each set contains indices of
    connected rings for that ring.
    """
    n_rings = len(rings)
    ring_sets = [set(r) for r in rings]
    spiro_set = set(spiro_atoms)
    
    connections = [set() for _ in range(n_rings)]
    
    for i in range(n_rings):
        for j in range(i + 1, n_rings):
            # Check if rings i and j share a spiro atom
            shared = ring_sets[i] & ring_sets[j] & spiro_set
            if shared:
                connections[i].add(j)
                connections[j].add(i)
    
    return connections


# =============================================================================
# P-24.5: INDICATED HYDROGEN IN SPIRO SYSTEMS
# =============================================================================

def get_indicated_hydrogen_spiro(numbering: Dict[int, int], mol) -> str:
    """
    Determine indicated hydrogen prefix for spiro system.
    
    Format: 2H,4''H-spiro[...]
    """
    from rdkit import Chem
    
    indicated_h = []
    
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in numbering:
            continue
        
        # Check if this is a position that needs indicated hydrogen
        # (typically sp3 nitrogen or certain ring positions)
        if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() > 0:
            pos = numbering[idx]
            indicated_h.append(pos)
    
    if indicated_h:
        return ','.join(f"{p}H" for p in sorted(indicated_h)) + '-'
    
    return ""
