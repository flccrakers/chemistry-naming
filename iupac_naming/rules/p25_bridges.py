"""
P-25.4: Bridged Fused Ring Systems
==================================

Rules for naming fused ring systems with additional bridges.

Reference: IUPAC Blue Book 2013, P-25.4

Bridged fused systems have:
1. A polycyclic base system
2. Additional bridge(s) connecting non-adjacent positions
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass


@dataclass  
class Bridge:
    """Represents a bridge in a fused system."""
    bridge_type: str          # 'epoxy', 'methano', 'ethano', etc.
    atom_indices: List[int]   # Atoms forming the bridge
    locant1: int              # First attachment position
    locant2: int              # Second attachment position
    heteroatoms: Dict[int, str]  # Position -> heteroatom in bridge


# =============================================================================
# BRIDGE TYPE IDENTIFICATION
# =============================================================================

BRIDGE_NAMES = {
    # Single atom bridges
    'O': 'epoxy',           # -O- bridge
    'N': 'epimino',         # -NH- bridge  
    'S': 'epithio',         # -S- bridge
    'CH2': 'methano',       # -CH2- bridge
    
    # Two atom bridges
    'CH2CH2': 'ethano',     # -CH2-CH2- bridge
    'CH2O': 'methanooxy',   # -CH2-O- bridge (methylene + oxygen)
    'OCH2': 'oxymethano',   # -O-CH2- bridge
    'OCH2O': 'methylenedioxy', # -O-CH2-O- bridge
    'NH': 'imino',          # =NH bridge (double bonded)
    'NCH3': 'methylimino',  # -N(CH3)- bridge
    
    # Three atom bridges
    'CH2CH2CH2': 'propano', # -CH2-CH2-CH2- bridge
    'OCH2CH2': 'oxoethano', # -O-CH2-CH2- bridge
    
    # Named bridges
    'benzo': 'benzo',       # Benzene ring as bridge
    'naphtho': 'naphtho',   # Naphthalene as bridge
}


def identify_bridge_type(atoms: List[int], mol) -> str:
    """
    Identify the type of bridge from its atoms.
    
    Args:
        atoms: Atom indices forming the bridge (excluding attachment points)
        mol: RDKit molecule
        
    Returns:
        Bridge type name
    """
    if not atoms:
        return "direct"  # Direct fusion (no bridging atoms)
    
    # Build signature of bridge atoms
    symbols = []
    for idx in atoms:
        atom = mol.GetAtomWithIdx(idx)
        sym = atom.GetSymbol()
        if sym == 'C':
            # Check for CH2, CH, C
            n_h = atom.GetTotalNumHs()
            if n_h >= 2:
                sym = 'CH2'
            elif n_h == 1:
                sym = 'CH'
        symbols.append(sym)
    
    signature = ''.join(symbols)
    
    # Check against known bridge types
    if signature in BRIDGE_NAMES:
        return BRIDGE_NAMES[signature]
    
    # Single atom bridges
    if len(atoms) == 1:
        sym = mol.GetAtomWithIdx(atoms[0]).GetSymbol()
        return BRIDGE_NAMES.get(sym, f"{sym.lower()}a")
    
    # Multi-atom bridges - try to construct name
    return construct_bridge_name(symbols)


def construct_bridge_name(symbols: List[str]) -> str:
    """
    Construct a bridge name from atom symbols.
    
    For chains: methano, ethano, propano, butano...
    With heteroatoms: oxa, aza, thia prefixes
    """
    n_carbons = sum(1 for s in symbols if s.startswith('C'))
    
    carbon_names = {
        1: 'methano',
        2: 'ethano', 
        3: 'propano',
        4: 'butano',
        5: 'pentano',
    }
    
    base = carbon_names.get(n_carbons, f'C{n_carbons}ano')
    
    # Add heteroatom prefixes
    het_prefixes = []
    het_map = {'O': 'oxa', 'N': 'aza', 'S': 'thia'}
    
    for i, sym in enumerate(symbols):
        if sym in het_map:
            het_prefixes.append(f"{i+1}-{het_map[sym]}")
    
    if het_prefixes:
        return ''.join(het_prefixes) + base
    
    return base


# =============================================================================
# BRIDGE LOCANT FORMATTING
# =============================================================================

def format_bridge_locants(locant1: int, locant2: int, 
                         prime1: str = "", prime2: str = "") -> str:
    """
    Format bridge attachment locants.
    
    Format: lower,higher where lower < higher numerically
    Primes are used for fusion locants from attached components.
    """
    # Ensure lower number first
    if locant1 <= locant2:
        return f"{locant1}{prime1},{locant2}{prime2}"
    else:
        return f"{locant2}{prime2},{locant1}{prime1}"


def format_bridge_prefix(bridge: Bridge) -> str:
    """
    Format complete bridge prefix.
    
    Format: locants-bridgetype
    Example: 9,12-epoxy
    """
    locs = format_bridge_locants(bridge.locant1, bridge.locant2)
    return f"{locs}-{bridge.bridge_type}"


# =============================================================================
# FINDING BRIDGES IN FUSED SYSTEMS
# =============================================================================

def find_bridges(mol, fused_system_atoms: Set[int]) -> List[Bridge]:
    """
    Find bridges in a fused polycyclic system.
    
    A bridge is identified when:
    1. Two atoms in the fused system are connected
    2. But they are not adjacent in the peripheral numbering
    3. The connection may be direct or through additional atoms
    
    Args:
        mol: RDKit molecule
        fused_system_atoms: Set of atom indices in the main fused system
        
    Returns:
        List of Bridge objects
    """
    from rdkit import Chem
    
    bridges = []
    ring_info = mol.GetRingInfo()
    
    # Find atoms that could be bridge endpoints
    # (atoms in fused system with connections outside normal ring fusion)
    
    for atom_idx in fused_system_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            
            # Check if this connection forms a bridge
            if n_idx not in fused_system_atoms:
                # This atom connects outside the fused system
                # Trace where it leads
                bridge_atoms, end_point = trace_bridge(
                    mol, n_idx, fused_system_atoms, {atom_idx}
                )
                
                if end_point is not None and end_point in fused_system_atoms:
                    # Found a bridge!
                    bridge_type = identify_bridge_type(bridge_atoms, mol)
                    
                    # Avoid duplicates (same bridge found from both ends)
                    if atom_idx < end_point:
                        bridges.append(Bridge(
                            bridge_type=bridge_type,
                            atom_indices=bridge_atoms,
                            locant1=atom_idx,
                            locant2=end_point,
                            heteroatoms={}
                        ))
    
    return bridges


def trace_bridge(mol, start_idx: int, fused_atoms: Set[int], 
                 exclude: Set[int]) -> Tuple[List[int], Optional[int]]:
    """
    Trace a potential bridge path.
    
    Args:
        mol: RDKit molecule
        start_idx: Starting atom (first bridge atom)
        fused_atoms: Set of fused system atoms
        exclude: Atoms already visited
        
    Returns:
        Tuple of (bridge_atoms, end_point) or ([], None) if not a bridge
    """
    visited = exclude.copy()
    path = []
    current = start_idx
    
    while current is not None:
        if current in visited:
            return ([], None)  # Cycle without reaching fused system
        
        visited.add(current)
        
        if current in fused_atoms and current not in exclude:
            # Reached the other end
            return (path, current)
        
        path.append(current)
        
        # Find next atom
        atom = mol.GetAtomWithIdx(current)
        next_atom = None
        
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                next_atom = n_idx
                break
        
        current = next_atom
    
    return ([], None)


# =============================================================================
# COMPLEX FUSION WITH BRIDGES
# =============================================================================

def name_bridged_fused_system(base_name: str, bridges: List[Bridge],
                              peripheral_numbering: Dict[int, int]) -> str:
    """
    Build complete name for bridged fused system.
    
    Format: bridge-prefixes + base-fusion-name
    Example: 9,12-epoxypyrido[...]phenazine
    
    Args:
        base_name: Name of the underlying fused system
        bridges: List of Bridge objects
        peripheral_numbering: Mapping of atom index to peripheral number
        
    Returns:
        Complete bridged fused system name
    """
    if not bridges:
        return base_name
    
    # Format each bridge prefix
    bridge_prefixes = []
    
    for bridge in bridges:
        # Convert atom indices to peripheral numbers
        loc1 = peripheral_numbering.get(bridge.locant1, '?')
        loc2 = peripheral_numbering.get(bridge.locant2, '?')
        
        # Sort locants
        if isinstance(loc1, int) and isinstance(loc2, int):
            if loc1 > loc2:
                loc1, loc2 = loc2, loc1
        
        prefix = f"{loc1},{loc2}-{bridge.bridge_type}"
        bridge_prefixes.append(prefix)
    
    # Sort bridge prefixes alphabetically by bridge type
    bridge_prefixes.sort(key=lambda x: x.split('-')[-1])
    
    return ''.join(bridge_prefixes) + base_name


# =============================================================================
# MULTI-LEVEL FUSION DESCRIPTORS
# =============================================================================

def format_fusion_descriptor(outer_locants: str, middle_locants: str,
                            inner_locants: str) -> str:
    """
    Format multi-level fusion descriptor for complex systems.
    
    Format: [outer:middle:inner]
    
    Example: [1'',2'':1',2'] means:
    - Outer component (double-prime) fused at positions 1,2
    - To middle component (single-prime) at positions 1,2
    
    Args:
        outer_locants: Locants with most primes
        middle_locants: Locants with fewer primes
        inner_locants: Locants closest to base
        
    Returns:
        Formatted fusion descriptor
    """
    parts = []
    
    if outer_locants:
        parts.append(outer_locants)
    if middle_locants:
        parts.append(middle_locants)
    if inner_locants:
        parts.append(inner_locants)
    
    return ':'.join(parts)


def assign_prime_levels(components: List[Dict], base_component: Dict) -> Dict[int, int]:
    """
    Assign prime levels to components based on distance from base.
    
    Base component: no primes
    First attached: single prime (')
    Second level: double prime ('')
    etc.
    
    Args:
        components: List of component dictionaries with 'atoms' sets
        base_component: The base component dictionary
        
    Returns:
        Dict mapping component index to prime level
    """
    prime_levels = {}
    
    # Base has level 0
    base_idx = components.index(base_component) if base_component in components else 0
    prime_levels[base_idx] = 0
    
    # BFS from base to assign levels
    visited = {base_idx}
    queue = [(base_idx, 0)]
    
    while queue:
        current_idx, level = queue.pop(0)
        current_atoms = components[current_idx]['atoms']
        
        for i, comp in enumerate(components):
            if i in visited:
                continue
            
            # Check if comp shares atoms with current
            if comp['atoms'] & current_atoms:
                prime_levels[i] = level + 1
                visited.add(i)
                queue.append((i, level + 1))
    
    return prime_levels


def format_prime_locant(locant: int, prime_level: int) -> str:
    """Format a locant with appropriate number of primes."""
    primes = "'" * prime_level
    return f"{locant}{primes}"
