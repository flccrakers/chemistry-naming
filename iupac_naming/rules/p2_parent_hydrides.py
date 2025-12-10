"""
P-2 PARENT HYDRIDES
===================

This module implements the rules for parent hydrides from Chapter P-2
of the IUPAC 2013 Blue Book.

Sections:
- P-21: Mononuclear and Acyclic Polynuclear Parent Hydrides
- P-22: Monocyclic Parent Hydrides (includes Hantzsch-Widman)
- P-23: Polyalicyclic Parent Hydrides (von Baeyer System)
- P-24: Spiro Ring Systems
- P-25: Fused and Bridged Fused Ring Systems
- P-26: Phane Nomenclature
- P-27: Fullerenes
- P-28: Ring Assemblies
- P-29: Substituent Group Prefixes
"""

from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass
from enum import Enum


# =============================================================================
# P-21 MONONUCLEAR AND ACYCLIC PARENT HYDRIDES
# =============================================================================

# P-21.1.1 Mononuclear parent hydrides with standard bonding numbers
MONONUCLEAR_HYDRIDES = {
    'B': 'borane',
    'C': 'methane',
    'Si': 'silane',
    'Ge': 'germane',
    'Sn': 'stannane',
    'Pb': 'plumbane',
    'N': 'azane',      # Also: ammonia (retained)
    'P': 'phosphane',  # Also: phosphine (retained)
    'As': 'arsane',
    'Sb': 'stibane',
    'Bi': 'bismuthane',
    'O': 'oxidane',    # Also: water (retained)
    'S': 'sulfane',
    'Se': 'selane',
    'Te': 'tellane',
}

# Retained names for common mononuclear hydrides (P-21.1.1)
RETAINED_MONONUCLEAR = {
    'NH3': 'ammonia',
    'PH3': 'phosphine',
    'H2O': 'water',
    'H2S': 'hydrogen sulfide',
}


def get_acyclic_hydride_name(element: str, n_atoms: int) -> str:
    """
    Get name for acyclic homogeneous parent hydride (P-21.2.2).
    
    Examples:
        C, 2 -> ethane
        C, 3 -> propane
        Si, 4 -> tetrasilane
    """
    # Hydrocarbon chain names (P-21.2.1)
    alkane_names = {
        1: 'methane', 2: 'ethane', 3: 'propane', 4: 'butane',
        5: 'pentane', 6: 'hexane', 7: 'heptane', 8: 'octane',
        9: 'nonane', 10: 'decane', 11: 'undecane', 12: 'dodecane',
        13: 'tridecane', 14: 'tetradecane', 15: 'pentadecane',
        16: 'hexadecane', 17: 'heptadecane', 18: 'octadecane',
        19: 'nonadecane', 20: 'icosane', 21: 'henicosane',
        22: 'docosane', 23: 'tricosane', 24: 'tetracosane',
        25: 'pentacosane', 26: 'hexacosane', 27: 'heptacosane',
        28: 'octacosane', 29: 'nonacosane', 30: 'triacontane'
    }
    
    if element == 'C':
        return alkane_names.get(n_atoms, f"n-{n_atoms}-ane")
    
    # For other elements, use prefix + element hydride name
    base = MONONUCLEAR_HYDRIDES.get(element, element.lower() + 'ane')
    if n_atoms == 1:
        return base
    
    from .p1_general import get_multiplicative_prefix
    prefix = get_multiplicative_prefix(n_atoms)
    
    # Remove 'ane' and add prefix + 'ane'
    if base.endswith('ane'):
        return prefix + base
    return prefix + base


# =============================================================================
# P-22 MONOCYCLIC PARENT HYDRIDES
# =============================================================================

# P-22.1.1 Saturated monocyclic hydrocarbons
def get_cycloalkane_name(n_atoms: int) -> str:
    """
    Get name for cycloalkane (P-22.1.1).
    """
    base = get_acyclic_hydride_name('C', n_atoms)
    return f"cyclo{base}"


# P-22.1.3 Retained names for monocyclic hydrocarbons
RETAINED_MONOCYCLIC_HYDROCARBONS = {
    6: 'benzene',       # Aromatic 6-membered
    5: 'cyclopentadiene',  # With unsaturation
}


# =============================================================================
# P-22.2 HANTZSCH-WIDMAN NOMENCLATURE
# =============================================================================

# P-22.2.2 'a' prefixes for heteroatoms (in order of priority)
HANTZSCH_WIDMAN_PREFIXES = {
    'F': 'fluora',
    'Cl': 'chlora',
    'Br': 'broma',
    'I': 'ioda',
    'O': 'oxa',
    'S': 'thia',
    'Se': 'selena',
    'Te': 'tellura',
    'N': 'aza',
    'P': 'phospha',
    'As': 'arsa',
    'Sb': 'stiba',
    'Bi': 'bisma',
    'Si': 'sila',
    'Ge': 'germa',
    'Sn': 'stanna',
    'Pb': 'plumba',
    'B': 'bora',
}

# P-22.2.2 Ring size stems for Hantzsch-Widman
HW_STEMS_UNSATURATED = {
    3: 'irene',   # or 'irine' with N
    4: 'ete',
    5: 'ole',
    6: 'ine',
    7: 'epine',
    8: 'ocine',
    9: 'onine',
    10: 'ecine',
}

HW_STEMS_SATURATED = {
    3: 'irane',    # or 'iridine' with N
    4: 'etane',    # or 'etidine' with N
    5: 'olane',    # or 'olidine' with N
    6: 'ane',      # or 'inane'/'idine' with N
    7: 'epane',
    8: 'ocane',
    9: 'onane',
    10: 'ecane',
}

# Special stems when nitrogen is present
HW_STEMS_WITH_NITROGEN_SAT = {
    3: 'iridine',
    4: 'etidine',
    5: 'olidine',
    6: 'inane',  # or perhydro- prefix with unsaturated name
}


def get_hantzsch_widman_name(ring_size: int, 
                             heteroatoms: Dict[int, str],
                             is_saturated: bool = True) -> str:
    """
    Generate Hantzsch-Widman name for a heterocycle (P-22.2.2).
    
    Args:
        ring_size: Number of atoms in the ring
        heteroatoms: Dict mapping position -> element symbol
        is_saturated: Whether the ring is saturated
    
    Returns:
        Hantzsch-Widman name
    """
    if ring_size < 3 or ring_size > 10:
        return None
    
    # Sort heteroatoms by priority (order in HANTZSCH_WIDMAN_PREFIXES)
    priority_order = list(HANTZSCH_WIDMAN_PREFIXES.keys())
    
    sorted_heteros = sorted(
        heteroatoms.items(),
        key=lambda x: priority_order.index(x[1]) if x[1] in priority_order else 999
    )
    
    # Build prefix
    from collections import Counter
    hetero_counts = Counter(e for _, e in sorted_heteros)
    
    from .p1_general import get_multiplicative_prefix
    
    prefix_parts = []
    for element in priority_order:
        if element in hetero_counts:
            count = hetero_counts[element]
            hw_prefix = HANTZSCH_WIDMAN_PREFIXES[element]
            mult_prefix = get_multiplicative_prefix(count)
            prefix_parts.append(mult_prefix + hw_prefix)
    
    prefix = ''.join(prefix_parts)
    
    # Get stem based on saturation and nitrogen presence
    has_nitrogen = 'N' in hetero_counts
    
    if is_saturated:
        if has_nitrogen and ring_size in HW_STEMS_WITH_NITROGEN_SAT:
            stem = HW_STEMS_WITH_NITROGEN_SAT[ring_size]
        else:
            stem = HW_STEMS_SATURATED.get(ring_size, 'ane')
    else:
        stem = HW_STEMS_UNSATURATED.get(ring_size, 'ine')
    
    # Elide vowel if needed
    from .p1_general import NameWriter
    name = NameWriter.elide_vowel(prefix, stem)
    
    return name


# =============================================================================
# P-22.2.1 RETAINED NAMES FOR HETEROCYCLES
# =============================================================================

# Monocyclic heterocycles with retained names
RETAINED_HETEROCYCLES = {
    # 5-membered with one heteroatom
    'furan': {'size': 5, 'heteroatoms': {'O': 1}, 'aromatic': True},
    'thiophene': {'size': 5, 'heteroatoms': {'S': 1}, 'aromatic': True},
    'pyrrole': {'size': 5, 'heteroatoms': {'N': 1}, 'aromatic': True},
    'selenophene': {'size': 5, 'heteroatoms': {'Se': 1}, 'aromatic': True},
    'tellurophene': {'size': 5, 'heteroatoms': {'Te': 1}, 'aromatic': True},
    
    # 5-membered with two heteroatoms
    'imidazole': {'size': 5, 'heteroatoms': {'N': 2}, 'aromatic': True, 'pattern': '1,3'},
    'pyrazole': {'size': 5, 'heteroatoms': {'N': 2}, 'aromatic': True, 'pattern': '1,2'},
    'oxazole': {'size': 5, 'heteroatoms': {'O': 1, 'N': 1}, 'aromatic': True, 'pattern': '1,3'},
    'isoxazole': {'size': 5, 'heteroatoms': {'O': 1, 'N': 1}, 'aromatic': True, 'pattern': '1,2'},
    'thiazole': {'size': 5, 'heteroatoms': {'S': 1, 'N': 1}, 'aromatic': True, 'pattern': '1,3'},
    'isothiazole': {'size': 5, 'heteroatoms': {'S': 1, 'N': 1}, 'aromatic': True, 'pattern': '1,2'},
    
    # 6-membered with one heteroatom
    'pyridine': {'size': 6, 'heteroatoms': {'N': 1}, 'aromatic': True},
    'pyran': {'size': 6, 'heteroatoms': {'O': 1}, 'aromatic': False},
    
    # 6-membered with two heteroatoms
    'pyrazine': {'size': 6, 'heteroatoms': {'N': 2}, 'aromatic': True, 'pattern': '1,4'},
    'pyrimidine': {'size': 6, 'heteroatoms': {'N': 2}, 'aromatic': True, 'pattern': '1,3'},
    'pyridazine': {'size': 6, 'heteroatoms': {'N': 2}, 'aromatic': True, 'pattern': '1,2'},
    
    # 6-membered with three heteroatoms
    'triazine': {'size': 6, 'heteroatoms': {'N': 3}, 'aromatic': True},
}


# =============================================================================
# P-25 FUSED RING SYSTEMS
# =============================================================================

# P-25.1.1 Retained names for fused hydrocarbon systems
RETAINED_FUSED_HYDROCARBONS = {
    'naphthalene': {'rings': 2, 'atoms': 10, 'pattern': '6,6'},
    'azulene': {'rings': 2, 'atoms': 10, 'pattern': '7,5'},
    'anthracene': {'rings': 3, 'atoms': 14, 'pattern': '6,6,6-linear'},
    'phenanthrene': {'rings': 3, 'atoms': 14, 'pattern': '6,6,6-angular'},
    'pyrene': {'rings': 4, 'atoms': 16},
    'fluorene': {'rings': 3, 'atoms': 13},
    'acenaphthylene': {'rings': 3, 'atoms': 12},
    'fluoranthene': {'rings': 4, 'atoms': 16},
    'triphenylene': {'rings': 4, 'atoms': 18},
    'chrysene': {'rings': 4, 'atoms': 18},
    'coronene': {'rings': 7, 'atoms': 24},
}

# P-25.2.1 Retained names for fused heterocyclic systems
RETAINED_FUSED_HETEROCYCLES = {
    # Benzofused 5-membered
    'indole': {'rings': 2, 'heteroatoms': {'N': 1}, 'base': 'benzene+pyrrole'},
    'benzofuran': {'rings': 2, 'heteroatoms': {'O': 1}, 'base': 'benzene+furan'},
    'benzothiophene': {'rings': 2, 'heteroatoms': {'S': 1}, 'base': 'benzene+thiophene'},
    'benzimidazole': {'rings': 2, 'heteroatoms': {'N': 2}, 'base': 'benzene+imidazole'},
    'benzoxazole': {'rings': 2, 'heteroatoms': {'O': 1, 'N': 1}, 'base': 'benzene+oxazole'},
    'benzothiazole': {'rings': 2, 'heteroatoms': {'S': 1, 'N': 1}, 'base': 'benzene+thiazole'},
    
    # Benzofused 6-membered
    'quinoline': {'rings': 2, 'heteroatoms': {'N': 1}, 'base': 'benzene+pyridine'},
    'isoquinoline': {'rings': 2, 'heteroatoms': {'N': 1}, 'base': 'benzene+pyridine'},
    'quinazoline': {'rings': 2, 'heteroatoms': {'N': 2}, 'base': 'benzene+pyrimidine'},
    'quinoxaline': {'rings': 2, 'heteroatoms': {'N': 2}, 'base': 'benzene+pyrazine'},
    'phthalazine': {'rings': 2, 'heteroatoms': {'N': 2}, 'base': 'benzene+pyridazine'},
    'cinnoline': {'rings': 2, 'heteroatoms': {'N': 2}},
    
    # Larger fused systems
    'acridine': {'rings': 3, 'heteroatoms': {'N': 1}},
    'phenazine': {'rings': 3, 'heteroatoms': {'N': 2}},
    'phenoxazine': {'rings': 3, 'heteroatoms': {'O': 1, 'N': 1}},
    'phenothiazine': {'rings': 3, 'heteroatoms': {'S': 1, 'N': 1}},
    
    # Purine and related
    'purine': {'rings': 2, 'heteroatoms': {'N': 4}, 'base': 'imidazole+pyrimidine'},
    'pteridine': {'rings': 2, 'heteroatoms': {'N': 4}},
    
    # Carbazole and xanthene
    'carbazole': {'rings': 3, 'heteroatoms': {'N': 1}},
    'xanthene': {'rings': 3, 'heteroatoms': {'O': 1}},
    'thioxanthene': {'rings': 3, 'heteroatoms': {'S': 1}},
}


# P-25.3 Fusion prefixes
FUSION_PREFIXES = {
    # From monocycles
    'benzene': 'benzo',
    'pyridine': 'pyrido',
    'pyrazine': 'pyrazino',
    'pyrimidine': 'pyrimido',
    'pyridazine': 'pyridazino',
    'furan': 'furo',
    'thiophene': 'thieno',
    'pyrrole': 'pyrrolo',
    'imidazole': 'imidazo',
    'pyrazole': 'pyrazolo',
    'oxazole': 'oxazolo',
    'thiazole': 'thiazolo',
    'isoxazole': 'isoxazolo',
    'isothiazole': 'isothiazolo',
    'triazole': 'triazolo',
    'tetrazole': 'tetrazolo',
    
    # From fused systems
    'naphthalene': 'naphtho',
    'quinoline': 'quinolino',
    'isoquinoline': 'isoquinolino',
    'indole': 'indolo',
    'purine': 'purino',
    'phenazine': 'phenazino',
}


# =============================================================================
# P-25.4 BRIDGES
# =============================================================================

# P-25.4.2.1 Bivalent bridge names
BIVALENT_BRIDGES = {
    # Hydrocarbon bridges
    'methano': {'atoms': 1, 'element': 'C'},
    'ethano': {'atoms': 2, 'element': 'C'},
    'propano': {'atoms': 3, 'element': 'C'},
    'butano': {'atoms': 4, 'element': 'C'},
    
    # Heteroatom bridges
    'epoxy': {'atoms': 1, 'element': 'O'},
    'epithio': {'atoms': 1, 'element': 'S'},
    'epimino': {'atoms': 1, 'element': 'N'},
    'epoxyethano': {'atoms': 2, 'elements': ['C', 'O']},
}


# =============================================================================
# P-28 RING ASSEMBLIES
# =============================================================================

# P-28.2 Prefixes for ring assemblies
RING_ASSEMBLY_PREFIXES = {
    2: 'bi',
    3: 'ter',
    4: 'quater',
    5: 'quinque',
    6: 'sexi',
    7: 'septi',
    8: 'octi',
    9: 'novi',
    10: 'deci',
}


def get_ring_assembly_name(base_name: str, n_units: int) -> str:
    """
    Get name for a ring assembly (P-28).
    
    Examples:
        biphenyl, terpyridine, quaterthiophene
    """
    if n_units < 2:
        return base_name
    
    prefix = RING_ASSEMBLY_PREFIXES.get(n_units, f"{n_units}-")
    
    # Handle special cases
    if base_name == 'phenyl' or base_name == 'benzene':
        return f"{prefix}phenyl"
    
    return f"{prefix}{base_name}"
