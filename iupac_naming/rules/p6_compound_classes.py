"""
P-6: Applications to Specific Classes of Compounds
===================================================

Nomenclature rules for specific compound classes.

Reference: IUPAC Blue Book 2013, Chapter P-6
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass


# =============================================================================
# P-61: HYDROCARBONS
# =============================================================================

@dataclass
class HydrocarbonType:
    """Classification of hydrocarbon types."""
    is_saturated: bool
    is_cyclic: bool
    is_aromatic: bool
    n_rings: int
    ring_sizes: List[int]


def classify_hydrocarbon(mol) -> HydrocarbonType:
    """
    Classify a hydrocarbon molecule.
    
    Reference: P-61
    """
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    is_cyclic = len(rings) > 0
    ring_sizes = [len(r) for r in rings]
    
    # Check saturation
    is_saturated = True
    is_aromatic = False
    
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() > 1:
            is_saturated = False
        if bond.GetIsAromatic():
            is_aromatic = True
    
    return HydrocarbonType(
        is_saturated=is_saturated,
        is_cyclic=is_cyclic,
        is_aromatic=is_aromatic,
        n_rings=len(rings),
        ring_sizes=ring_sizes
    )


# =============================================================================
# P-62: AMINES AND IMINES
# =============================================================================

AMINE_RETAINED_NAMES = {
    'methylamine': 'CN',
    'ethylamine': 'CCN',
    'propylamine': 'CCCN',
    'aniline': 'Nc1ccccc1',
    'benzylamine': 'NCc1ccccc1',
}


def name_amine(mol, primary_n_idx: int) -> str:
    """
    Generate name for an amine.
    
    Reference: P-62
    """
    from rdkit import Chem
    
    n_atom = mol.GetAtomWithIdx(primary_n_idx)
    
    # Count hydrogens to determine amine type
    n_h = n_atom.GetTotalNumHs()
    
    if n_h == 2:
        return "primary amine"
    elif n_h == 1:
        return "secondary amine"
    else:
        return "tertiary amine"


# =============================================================================
# P-63: ALCOHOLS AND PHENOLS
# =============================================================================

def name_alcohol(mol, oh_oxygen_idx: int, chain_length: int) -> Tuple[str, str]:
    """
    Generate alcohol name components.
    
    Reference: P-63
    
    Returns:
        Tuple of (prefix, suffix)
    """
    # Determine position of OH
    o_atom = mol.GetAtomWithIdx(oh_oxygen_idx)
    
    # Find the carbon attached to OH
    c_neighbor = None
    for neighbor in o_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
            c_neighbor = neighbor.GetIdx()
            break
    
    prefix = ""
    suffix = "ol"
    
    return (prefix, suffix)


# =============================================================================
# P-64: KETONES
# =============================================================================

def name_ketone(mol, carbonyl_c_idx: int) -> Tuple[str, str]:
    """
    Generate ketone name components.
    
    Reference: P-64
    """
    prefix = ""
    suffix = "one"
    
    return (prefix, suffix)


KETONE_RETAINED_NAMES = {
    'acetone': 'CC(=O)C',
    'benzophenone': 'c1ccc(cc1)C(=O)c2ccccc2',
}


# =============================================================================
# P-65: CARBOXYLIC ACIDS
# =============================================================================

ACID_RETAINED_NAMES = {
    'formic acid': 'C(=O)O',
    'acetic acid': 'CC(=O)O',
    'propionic acid': 'CCC(=O)O',
    'butyric acid': 'CCCC(=O)O',
    'valeric acid': 'CCCCC(=O)O',
    'benzoic acid': 'c1ccc(cc1)C(=O)O',
    'oxalic acid': 'C(=O)(C(=O)O)O',
}


def name_carboxylic_acid(mol, cooh_c_idx: int) -> str:
    """
    Generate carboxylic acid name.
    
    Reference: P-65
    """
    return "oic acid"


# =============================================================================
# P-66: ALDEHYDES AND THEIR DERIVATIVES
# =============================================================================

ALDEHYDE_RETAINED_NAMES = {
    'formaldehyde': 'C=O',
    'acetaldehyde': 'CC=O',
    'benzaldehyde': 'c1ccc(cc1)C=O',
}


def name_aldehyde(chain_length: int) -> str:
    """Generate aldehyde suffix."""
    return "al"


# =============================================================================
# P-67: NITRILES, ISOCYANIDES
# =============================================================================

def name_nitrile(chain_length: int) -> str:
    """Generate nitrile suffix."""
    return "nitrile"


def name_carbonitrile() -> str:
    """For -CN attached to ring."""
    return "carbonitrile"


# =============================================================================
# P-68: AMIDES, HYDRAZIDES, IMIDES
# =============================================================================

AMIDE_RETAINED_NAMES = {
    'formamide': 'C(=O)N',
    'acetamide': 'CC(=O)N',
    'benzamide': 'c1ccc(cc1)C(=O)N',
}


def name_amide(is_primary: bool) -> str:
    """Generate amide suffix."""
    if is_primary:
        return "amide"
    return "carboxamide"
