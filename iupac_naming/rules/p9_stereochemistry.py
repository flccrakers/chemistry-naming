"""
P-9: Stereochemistry
====================

Rules for stereochemical nomenclature including:
- E/Z configuration for double bonds
- R/S configuration for tetrahedral centers
- cis/trans for cyclic systems
- Other stereochemical descriptors

Reference: IUPAC Blue Book 2013, Chapter P-9
"""

from typing import List, Dict, Optional, Tuple, Set
from enum import Enum, auto


class StereoDescriptor(Enum):
    """Stereochemical descriptors."""
    R = auto()      # Rectus (clockwise)
    S = auto()      # Sinister (counterclockwise)
    E = auto()      # Entgegen (opposite)
    Z = auto()      # Zusammen (together)
    CIS = auto()    # Same side
    TRANS = auto()  # Opposite side
    ENDO = auto()   # Inside
    EXO = auto()    # Outside
    SYN = auto()    # Same side (general)
    ANTI = auto()   # Opposite side (general)
    ALPHA = auto()  # Below plane
    BETA = auto()   # Above plane


# =============================================================================
# P-91: CIP PRIORITY RULES
# =============================================================================

# Cahn-Ingold-Prelog priority order (higher atomic number = higher priority)
CIP_PRIORITY_BY_ATOM = {
    'I': 53,
    'Br': 35,
    'Cl': 17,
    'S': 16,
    'P': 15,
    'F': 9,
    'O': 8,
    'N': 7,
    'C': 6,
    'B': 5,
    'H': 1,
}


def get_cip_priority(atom_symbol: str) -> int:
    """Get CIP priority for an atom."""
    return CIP_PRIORITY_BY_ATOM.get(atom_symbol, 0)


def assign_cip_priorities(mol, center_idx: int) -> List[Tuple[int, int]]:
    """
    Assign CIP priorities to substituents around a stereocenter.
    
    Reference: P-91.1
    
    Args:
        mol: RDKit molecule
        center_idx: Index of stereocenter atom
        
    Returns:
        List of (neighbor_idx, priority) tuples, highest priority first
    """
    center = mol.GetAtomWithIdx(center_idx)
    neighbors = [(n.GetIdx(), get_cip_priority(n.GetSymbol())) 
                 for n in center.GetNeighbors()]
    
    # Sort by priority (descending)
    return sorted(neighbors, key=lambda x: x[1], reverse=True)


# =============================================================================
# P-92: R/S CONFIGURATION
# =============================================================================

def determine_rs_configuration(mol, chiral_center_idx: int) -> Optional[StereoDescriptor]:
    """
    Determine R/S configuration at a tetrahedral stereocenter.
    
    Reference: P-92
    
    Args:
        mol: RDKit molecule
        chiral_center_idx: Index of the chiral center
        
    Returns:
        StereoDescriptor.R or StereoDescriptor.S, or None if not determinable
    """
    from rdkit import Chem
    
    atom = mol.GetAtomWithIdx(chiral_center_idx)
    chiral_tag = atom.GetChiralTag()
    
    if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        return StereoDescriptor.R
    elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return StereoDescriptor.S
    
    return None


def format_rs_descriptor(descriptor: StereoDescriptor, locant: int) -> str:
    """Format R/S descriptor with locant."""
    if descriptor == StereoDescriptor.R:
        return f"({locant}R)"
    elif descriptor == StereoDescriptor.S:
        return f"({locant}S)"
    return ""


# =============================================================================
# P-93: E/Z CONFIGURATION
# =============================================================================

def determine_ez_configuration(mol, double_bond_idx: int) -> Optional[StereoDescriptor]:
    """
    Determine E/Z configuration at a double bond.
    
    Reference: P-93
    
    Args:
        mol: RDKit molecule
        double_bond_idx: Index of the double bond
        
    Returns:
        StereoDescriptor.E or StereoDescriptor.Z, or None
    """
    from rdkit import Chem
    
    bond = mol.GetBondWithIdx(double_bond_idx)
    stereo = bond.GetStereo()
    
    if stereo == Chem.BondStereo.STEREOE:
        return StereoDescriptor.E
    elif stereo == Chem.BondStereo.STEREOZ:
        return StereoDescriptor.Z
    
    return None


def format_ez_descriptor(descriptor: StereoDescriptor, locant: int) -> str:
    """Format E/Z descriptor with locant."""
    if descriptor == StereoDescriptor.E:
        return f"({locant}E)"
    elif descriptor == StereoDescriptor.Z:
        return f"({locant}Z)"
    return ""


# =============================================================================
# P-94: CIS/TRANS FOR CYCLIC SYSTEMS
# =============================================================================

def determine_cis_trans_ring(mol, ring_atoms: List[int], 
                             sub1_idx: int, sub2_idx: int) -> Optional[StereoDescriptor]:
    """
    Determine cis/trans relationship for substituents on a ring.
    
    Reference: P-94
    
    Args:
        mol: RDKit molecule
        ring_atoms: List of atom indices in the ring
        sub1_idx: First substituent attachment point
        sub2_idx: Second substituent attachment point
        
    Returns:
        StereoDescriptor.CIS or StereoDescriptor.TRANS, or None
    """
    # This requires 3D coordinates or explicit stereo bonds
    # Simplified implementation
    return None


# =============================================================================
# P-95: POLYCYCLIC STEREOCHEMISTRY
# =============================================================================

def determine_endo_exo(mol, bridge_atom_idx: int, 
                       reference_bridge: List[int]) -> Optional[StereoDescriptor]:
    """
    Determine endo/exo relationship in bridged systems.
    
    Reference: P-95
    """
    return None


# =============================================================================
# P-96: PRIMED LOCANTS FOR STEREOCHEMISTRY
# =============================================================================

def get_stereo_locants(descriptors: List[Tuple[int, StereoDescriptor]]) -> str:
    """
    Format multiple stereochemical descriptors.
    
    Reference: P-96
    
    Example: (2R,3S)-butane-2,3-diol
    """
    parts = []
    for locant, descriptor in sorted(descriptors):
        if descriptor == StereoDescriptor.R:
            parts.append(f"{locant}R")
        elif descriptor == StereoDescriptor.S:
            parts.append(f"{locant}S")
        elif descriptor == StereoDescriptor.E:
            parts.append(f"{locant}E")
        elif descriptor == StereoDescriptor.Z:
            parts.append(f"{locant}Z")
    
    if parts:
        return f"({','.join(parts)})-"
    return ""


# =============================================================================
# STEREOCHEMISTRY ANALYSIS
# =============================================================================

def find_stereocenters(mol) -> List[int]:
    """
    Find all stereocenters in a molecule.
    
    Returns:
        List of atom indices that are stereocenters
    """
    from rdkit import Chem
    
    stereocenters = []
    
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    for atom in mol.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            stereocenters.append(atom.GetIdx())
    
    return stereocenters


def find_stereogenic_double_bonds(mol) -> List[int]:
    """
    Find all stereogenic double bonds in a molecule.
    
    Returns:
        List of bond indices that are stereogenic double bonds
    """
    from rdkit import Chem
    
    stereobonds = []
    
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo != Chem.BondStereo.STEREONONE:
                stereobonds.append(bond.GetIdx())
    
    return stereobonds


def get_full_stereo_prefix(mol) -> str:
    """
    Generate complete stereochemical prefix for a molecule.
    
    Example: "(2R,3S,5E)-"
    """
    descriptors = []
    
    # R/S centers
    for center_idx in find_stereocenters(mol):
        config = determine_rs_configuration(mol, center_idx)
        if config:
            # Get locant (simplified - would need proper numbering)
            locant = center_idx + 1
            descriptors.append((locant, config))
    
    # E/Z bonds
    for bond_idx in find_stereogenic_double_bonds(mol):
        config = determine_ez_configuration(mol, bond_idx)
        if config:
            bond = mol.GetBondWithIdx(bond_idx)
            locant = min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) + 1
            descriptors.append((locant, config))
    
    return get_stereo_locants(descriptors)
