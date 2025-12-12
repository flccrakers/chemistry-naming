"""
P-5: Selecting Preferred IUPAC Names (PIN)
==========================================

Rules for selecting preferred IUPAC names when multiple valid names exist.

Reference: IUPAC Blue Book 2013, Chapter P-5
"""

from typing import List, Dict, Optional, Tuple
from enum import Enum, auto


class NameType(Enum):
    """Types of IUPAC names."""
    PIN = auto()           # Preferred IUPAC Name
    RETAINED = auto()      # Retained traditional name
    SYSTEMATIC = auto()    # Systematic name
    TRIVIAL = auto()       # Trivial/common name
    ACCEPTABLE = auto()    # Acceptable alternative


# =============================================================================
# P-51: GENERAL PRINCIPLES FOR NAME SELECTION
# =============================================================================

class NamePriority:
    """
    Priority rules for selecting between alternative names.
    
    Reference: P-51
    """
    
    # Priority order for selecting principal characteristic group
    # Higher index = higher priority
    PRINCIPAL_GROUP_PRIORITY = [
        'radicals',           # Lowest
        'anions',
        'cations', 
        'acids',              # carboxylic > sulfonic > sulfinic etc.
        'anhydrides',
        'esters',
        'acid_halides',
        'amides',
        'hydrazides',
        'imides',
        'nitriles',
        'aldehydes',
        'ketones',
        'alcohols',
        'hydroperoxides',
        'amines',
        'imines',
        'ethers',
        'peroxides',
    ]
    
    @classmethod
    def get_priority(cls, group_type: str) -> int:
        """Get priority index for a characteristic group."""
        try:
            return cls.PRINCIPAL_GROUP_PRIORITY.index(group_type)
        except ValueError:
            return -1


# =============================================================================
# P-52: SELECTING PREFERRED RETAINED NAMES
# =============================================================================

# Compounds that MUST use retained names as PIN
MANDATORY_RETAINED_NAMES = {
    # P-52.1.1 - Saturated monocyclic hydrocarbons
    'cyclopropane', 'cyclobutane', 'cyclopentane', 'cyclohexane',
    'cycloheptane', 'cyclooctane',
    
    # P-52.2.1 - Unsaturated hydrocarbons
    'benzene', 'naphthalene', 'anthracene', 'phenanthrene',
    'pyrene', 'chrysene', 'coronene', 'azulene',
    
    # P-52.2.2 - Heterocycles
    'furan', 'thiophene', 'pyrrole', 'pyridine',
    'pyrazine', 'pyrimidine', 'pyridazine', 'pyrazole',
    'imidazole', 'oxazole', 'thiazole', 'isoxazole', 'isothiazole',
    'quinoline', 'isoquinoline', 'indole', 'benzofuran', 'benzothiophene',
    'purine', 'pteridine', 'phenazine', 'carbazole', 'acridine',
    
    # P-52.2.3 - Functional compounds
    'formic acid', 'acetic acid', 'oxalic acid',
    'formaldehyde', 'acetaldehyde',
    'acetone', 'benzophenone',
    'methanol', 'ethanol', 'phenol',
    'aniline', 'glycine',
}


# =============================================================================
# P-53: SELECTING THE PRINCIPAL CHAIN
# =============================================================================

def select_principal_chain(chains: List[List[int]], mol) -> List[int]:
    """
    Select the principal chain according to IUPAC rules.
    
    Priority (P-53.1):
    1. Maximum number of substituents corresponding to principal characteristic group
    2. Maximum number of substituents corresponding to suffixes
    3. Maximum length
    4. Maximum number of heteroatoms (for chains)
    5. Maximum number of multiple bonds
    6. Maximum number of double bonds
    7. Lowest locants for principal characteristic groups
    8. Lowest locants for suffixes
    9. Maximum number of prefix substituents
    10. Lowest locants for prefix substituents
    
    Args:
        chains: List of atom index lists representing possible chains
        mol: RDKit molecule object
        
    Returns:
        Selected principal chain
    """
    if not chains:
        return []
    
    if len(chains) == 1:
        return chains[0]
    
    def chain_score(chain):
        """Calculate priority score for a chain."""
        length = len(chain)
        
        # Count heteroatoms in chain
        heteroatoms = sum(1 for idx in chain 
                        if mol.GetAtomWithIdx(idx).GetSymbol() != 'C')
        
        # Count multiple bonds
        multiple_bonds = 0
        double_bonds = 0
        for i in range(len(chain) - 1):
            bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
            if bond:
                bt = bond.GetBondTypeAsDouble()
                if bt > 1:
                    multiple_bonds += 1
                    if bt == 2:
                        double_bonds += 1
        
        # Return tuple for comparison (higher = better)
        return (length, heteroatoms, multiple_bonds, double_bonds)
    
    return max(chains, key=chain_score)


# =============================================================================
# P-54: SELECTING THE SENIOR RING SYSTEM
# =============================================================================

def select_senior_ring_system(ring_systems: List[Dict], mol) -> Dict:
    """
    Select the senior ring system according to IUPAC rules.
    
    Priority (P-54):
    1. Rings containing nitrogen
    2. Ring with most rings
    3. Ring with most atoms
    4. Ring with most heteroatoms
    5. Ring with most variety of heteroatoms
    6. Ring with larger individual rings
    
    Args:
        ring_systems: List of ring system dictionaries
        mol: RDKit molecule object
        
    Returns:
        Selected senior ring system
    """
    if not ring_systems:
        return {}
    
    if len(ring_systems) == 1:
        return ring_systems[0]
    
    def system_score(system):
        atoms = system.get('atoms', set())
        rings = system.get('rings', [])
        
        # Check for nitrogen
        has_nitrogen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'N' 
                          for idx in atoms)
        
        # Count heteroatoms
        heteroatoms = []
        for idx in atoms:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            if sym != 'C':
                heteroatoms.append(sym)
        
        n_heteroatoms = len(heteroatoms)
        variety = len(set(heteroatoms))
        
        # Max ring size
        max_ring_size = max(len(r) for r in rings) if rings else 0
        
        return (
            int(has_nitrogen),
            len(rings),
            len(atoms),
            n_heteroatoms,
            variety,
            max_ring_size
        )
    
    return max(ring_systems, key=system_score)


# =============================================================================
# P-55: NUMBERING
# =============================================================================

def apply_lowest_locant_rule(locants: List[int], alternative: List[int]) -> List[int]:
    """
    Apply the lowest locant rule (P-55.1.1).
    
    Compare locant sets and return the one that is "lower" at the
    first point of difference.
    
    Args:
        locants: First set of locants
        alternative: Alternative set of locants
        
    Returns:
        The lower set of locants
    """
    # Sort both lists
    loc1 = sorted(locants)
    loc2 = sorted(alternative)
    
    # Compare element by element
    for l1, l2 in zip(loc1, loc2):
        if l1 < l2:
            return locants
        elif l2 < l1:
            return alternative
    
    # If equal, prefer the first (original)
    return locants


def calculate_optimal_numbering(atoms: List[int], mol, 
                                characteristic_groups: List[int] = None,
                                unsaturation: List[int] = None) -> Dict[int, int]:
    """
    Calculate optimal numbering for a chain or ring.
    
    Reference: P-55
    
    Args:
        atoms: List of atom indices in order
        mol: RDKit molecule
        characteristic_groups: Atom indices with principal groups
        unsaturation: Atom indices with double/triple bonds
        
    Returns:
        Mapping of atom index to position number
    """
    n = len(atoms)
    
    if n == 0:
        return {}
    
    # Forward numbering (1, 2, 3, ...)
    forward = {atoms[i]: i + 1 for i in range(n)}
    
    # Reverse numbering (n, n-1, ..., 1)
    reverse = {atoms[i]: n - i for i in range(n)}
    
    # Determine which is preferred
    if characteristic_groups:
        # Get locants for both numberings
        fwd_locants = sorted([forward[idx] for idx in characteristic_groups if idx in forward])
        rev_locants = sorted([reverse[idx] for idx in characteristic_groups if idx in reverse])
        
        # Compare
        if fwd_locants <= rev_locants:
            return forward
        else:
            return reverse
    
    elif unsaturation:
        fwd_locants = sorted([forward[idx] for idx in unsaturation if idx in forward])
        rev_locants = sorted([reverse[idx] for idx in unsaturation if idx in reverse])
        
        if fwd_locants <= rev_locants:
            return forward
        else:
            return reverse
    
    return forward


# =============================================================================
# P-56: NAME CONSTRUCTION ORDER
# =============================================================================

NAME_COMPONENT_ORDER = [
    'detachable_prefixes',      # e.g., N-methyl
    'nondetachable_prefixes',   # e.g., cyclo
    'multiplicative_prefixes',  # e.g., di, tri
    'substituent_prefixes',     # Alphabetically ordered
    'parent_hydride',           # Main structure name
    'unsaturation_suffix',      # -ene, -yne
    'principal_suffix',         # -ol, -one, etc.
    'functional_suffix',        # -oic acid, etc.
]


def assemble_name(components: Dict[str, str]) -> str:
    """
    Assemble a complete IUPAC name from components.
    
    Reference: P-56
    
    Args:
        components: Dictionary of name components
        
    Returns:
        Complete IUPAC name
    """
    parts = []
    
    for component_type in NAME_COMPONENT_ORDER:
        if component_type in components and components[component_type]:
            parts.append(components[component_type])
    
    return ''.join(parts)
