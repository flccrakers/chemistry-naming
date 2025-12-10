"""
P-4 RULES FOR NAME CONSTRUCTION
===============================

This module implements the seniority and selection rules from Chapter P-4
of the IUPAC 2013 Blue Book.

Sections:
- P-41: Seniority Order for Classes
- P-42: Seniority Order for Acids
- P-43: Seniority Order for Suffixes
- P-44: Seniority Order for Parent Structures
- P-45: Selection of Preferred IUPAC Names
- P-46: Principal Chain in Substituent Groups
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass
from enum import Enum, auto
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_structures import (
    RingInfo, CharacteristicGroup, CharacteristicGroupClass,
    FusedRingSystem, SpiroSystem, VonBaeyerSystem
)


# =============================================================================
# P-41 SENIORITY ORDER FOR CLASSES
# =============================================================================

# Classes in decreasing order of seniority (P-41)
# Higher number = higher seniority
CLASS_SENIORITY = {
    # Class 1-4: Radicals, anions, cations, zwitterions
    'radical': 100,
    'anion': 99,
    'cation': 98,
    'zwitterion': 97,
    
    # Class 5-6: Acids
    'carboxylic_acid': 90,
    'carbothioic_acid': 89,
    'carboselenoic_acid': 88,
    'sulfonic_acid': 87,
    'sulfinic_acid': 86,
    'selenonic_acid': 85,
    'seleninic_acid': 84,
    'phosphonic_acid': 83,
    
    # Class 7: Anhydrides and esters
    'anhydride': 80,
    'ester': 79,
    'lactone': 78,
    
    # Class 8: Acid halides
    'acyl_halide': 75,
    
    # Class 9-11: Amides and related
    'amide': 70,
    'lactam': 69,
    'imide': 68,
    'hydrazide': 67,
    
    # Class 12-13: Nitriles and related
    'nitrile': 65,
    'isocyanide': 64,
    
    # Class 14: Aldehydes
    'aldehyde': 60,
    
    # Class 15: Ketones
    'ketone': 55,
    
    # Class 16: Alcohols and phenols
    'alcohol': 50,
    'phenol': 49,
    
    # Class 17: Thiols and analogs
    'thiol': 45,
    'selenol': 44,
    
    # Class 18: Hydroperoxides
    'hydroperoxide': 42,
    
    # Class 19: Amines
    'amine': 40,
    'imine': 39,
    
    # Class 20: Ethers and peroxides
    'ether': 35,
    'peroxide': 34,
    
    # Lowest: Parent hydrides
    'hydrocarbon': 10,
    'heterocycle': 10,
}


def get_class_seniority(class_name: str) -> int:
    """
    Get the seniority value for a compound class (P-41).
    
    Higher value = higher seniority = named as suffix
    """
    return CLASS_SENIORITY.get(class_name.lower(), 0)


def compare_class_seniority(class1: str, class2: str) -> int:
    """
    Compare two classes for seniority.
    
    Returns:
        -1 if class1 is senior (higher priority)
        0 if equal
        1 if class2 is senior
    """
    s1 = get_class_seniority(class1)
    s2 = get_class_seniority(class2)
    
    if s1 > s2:
        return -1
    elif s1 < s2:
        return 1
    return 0


# =============================================================================
# P-44 SENIORITY ORDER FOR PARENT STRUCTURES
# =============================================================================

# P-44.2 Seniority for rings and ring systems
def compare_ring_seniority(ring1: RingInfo, ring2: RingInfo) -> int:
    """
    Compare two rings for seniority (P-44.2).
    
    Criteria (in order):
    1. Heterocycle > carbocycle (P-44.2.1.2)
    2. Contains nitrogen > no nitrogen (P-44.2.1.3)  
    3. Contains O, S, Se, Te, N, P, As, etc. in that order (P-44.2.1.4)
    4. Greater number of rings (P-44.2.1.5)
    5. Greater number of skeletal atoms (P-44.2.1.6)
    6. Greater number of heteroatoms (P-44.2.1.7)
    
    Returns:
        -1 if ring1 is senior
        0 if equal
        1 if ring2 is senior
    """
    # Criterion 1: Heterocycle > carbocycle
    if ring1.is_heterocyclic and not ring2.is_heterocyclic:
        return -1
    if ring2.is_heterocyclic and not ring1.is_heterocyclic:
        return 1
    
    # Criterion 2: Contains nitrogen
    if ring1.has_nitrogen and not ring2.has_nitrogen:
        return -1
    if ring2.has_nitrogen and not ring1.has_nitrogen:
        return 1
    
    # Criterion 3: Heteroatom priority
    hetero_priority = ['O', 'S', 'Se', 'Te', 'N', 'P', 'As', 'Sb', 'Bi', 
                       'Si', 'Ge', 'Sn', 'Pb', 'B']
    
    def get_senior_heteroatom(heteroatoms: Dict[int, str]) -> int:
        if not heteroatoms:
            return 999
        elements = set(heteroatoms.values())
        for i, elem in enumerate(hetero_priority):
            if elem in elements:
                return i
        return 999
    
    h1 = get_senior_heteroatom(ring1.heteroatoms)
    h2 = get_senior_heteroatom(ring2.heteroatoms)
    
    if h1 < h2:
        return -1
    if h1 > h2:
        return 1
    
    # Criterion 5: Greater number of skeletal atoms
    if ring1.size > ring2.size:
        return -1
    if ring1.size < ring2.size:
        return 1
    
    # Criterion 6: Greater number of heteroatoms
    if ring1.hetero_count > ring2.hetero_count:
        return -1
    if ring1.hetero_count < ring2.hetero_count:
        return 1
    
    return 0


# P-44.2.2.2.3 Seniority for fused ring systems
def compare_fused_system_seniority(sys1: FusedRingSystem, 
                                   sys2: FusedRingSystem) -> int:
    """
    Compare two fused ring systems for seniority (P-44.2.2.2.3).
    
    Criteria:
    (a) Larger individual ring at first point of difference
    (b) Greater number of rings in horizontal row
    (c) Lower letters in fusion descriptor
    (d) Lower numbers in fusion descriptor
    (e) Senior ring component
    """
    # Compare by number of rings
    if len(sys1.rings) > len(sys2.rings):
        return -1
    if len(sys1.rings) < len(sys2.rings):
        return 1
    
    # Compare by total atoms
    if len(sys1.all_atoms) > len(sys2.all_atoms):
        return -1
    if len(sys1.all_atoms) < len(sys2.all_atoms):
        return 1
    
    # Compare individual ring sizes (sorted descending)
    sizes1 = sorted([r.size for r in sys1.rings], reverse=True)
    sizes2 = sorted([r.size for r in sys2.rings], reverse=True)
    
    for s1, s2 in zip(sizes1, sizes2):
        if s1 > s2:
            return -1
        if s1 < s2:
            return 1
    
    return 0


# =============================================================================
# P-44.3 SENIORITY OF ACYCLIC CHAINS (Principal Chain)
# =============================================================================

def select_principal_chain(chains: List[List[int]], 
                          atom_info: Dict[int, dict]) -> List[int]:
    """
    Select the principal chain from candidates (P-44.3).
    
    Criteria (in order):
    1. Maximum number of principal characteristic groups
    2. Maximum length
    3. Maximum number of multiple bonds
    4. Maximum number of double bonds
    5. Lowest locants for principal characteristic groups
    6. Lowest locants for multiple bonds
    7. Lowest locants for double bonds
    8. Maximum number of substituents
    """
    if not chains:
        return []
    
    if len(chains) == 1:
        return chains[0]
    
    # For now, select longest chain
    return max(chains, key=len)


# =============================================================================
# P-44.4 ADDITIONAL SENIORITY CRITERIA
# =============================================================================

def apply_low_locant_rule(locant_sets: List[List[int]]) -> int:
    """
    Apply the lowest locants rule (P-44.4).
    
    Returns index of the set with lowest locants.
    """
    from .p1_general import compare_locant_sets
    
    best_idx = 0
    for i in range(1, len(locant_sets)):
        if compare_locant_sets(locant_sets[i], locant_sets[best_idx]) < 0:
            best_idx = i
    
    return best_idx


# =============================================================================
# P-45 SELECTION OF PREFERRED IUPAC NAME
# =============================================================================

class PINSelectionCriteria(Enum):
    """
    Criteria for selecting Preferred IUPAC Names (P-45).
    """
    # P-45.1 Multiplication
    MULTIPLICATIVE_OVER_SUBSTITUTIVE = auto()
    
    # P-45.2 Substituent criteria
    MAX_SUBSTITUENTS = auto()
    LOWER_LOCANTS_FOR_SUBSTITUENTS = auto()
    
    # P-45.5 Alphanumerical order
    ALPHANUMERICAL = auto()
    
    # P-45.6 Stereochemistry
    STEREOCHEMISTRY = auto()


def select_preferred_name(candidates: List[str], 
                         criteria: List[PINSelectionCriteria] = None) -> str:
    """
    Select the Preferred IUPAC Name from candidates (P-45).
    """
    if not candidates:
        return ""
    
    if len(candidates) == 1:
        return candidates[0]
    
    # Apply default criteria
    if criteria is None:
        criteria = [
            PINSelectionCriteria.MULTIPLICATIVE_OVER_SUBSTITUTIVE,
            PINSelectionCriteria.MAX_SUBSTITUENTS,
            PINSelectionCriteria.LOWER_LOCANTS_FOR_SUBSTITUENTS,
            PINSelectionCriteria.ALPHANUMERICAL,
        ]
    
    # For now, return first candidate
    # Full implementation would apply all criteria
    return candidates[0]


# =============================================================================
# P-46 PRINCIPAL CHAIN IN SUBSTITUENT GROUPS
# =============================================================================

def select_principal_substituent_chain(chains: List[List[int]],
                                       atom_info: Dict[int, dict]) -> List[int]:
    """
    Select principal chain in a substituent group (P-46).
    
    Criteria (P-46.1):
    1. Greater number of heteroatoms
    2. Greater number of skeletal atoms
    3. Greater number of heteroatoms in set order
    4. Greater number of multiple bonds
    5. One or more atoms with nonstandard bonding numbers
    6. Lowest locants for heteroatoms
    7. Lowest locants for specific heteroatoms in order
    8. Lowest locants for free valences
    9. Lowest locants for multiple bonds
    10. Lowest locants for nonstandard bonding
    11. Greatest number of substituents
    12. Lowest locants for substituents
    13. Lowest locants for substituents cited first
    """
    if not chains:
        return []
    
    if len(chains) == 1:
        return chains[0]
    
    # Simplified: return longest chain
    return max(chains, key=len)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def determine_parent_structure(rings: List[RingInfo],
                               chains: List[List[int]],
                               characteristic_groups: List[CharacteristicGroup]) -> str:
    """
    Determine the type of parent structure for the molecule.
    
    Returns one of:
    - 'fused_ring'
    - 'spiro'
    - 'bridged'
    - 'ring_assembly'
    - 'ring_with_chain'
    - 'acyclic'
    """
    if not rings:
        return 'acyclic'
    
    if len(rings) == 1:
        if chains:
            return 'ring_with_chain'
        return 'monocyclic'
    
    # Check for fused rings (share 2+ atoms)
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            shared = set(rings[i].atoms) & set(rings[j].atoms)
            if len(shared) >= 2:
                return 'fused_ring'
            elif len(shared) == 1:
                return 'spiro'
    
    return 'ring_assembly'
