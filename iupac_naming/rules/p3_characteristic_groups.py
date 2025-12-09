"""
P-3 CHARACTERISTIC (FUNCTIONAL) AND SUBSTITUENT GROUPS
======================================================

This module implements rules for functional groups and substituents
from Chapter P-3 of the IUPAC 2013 Blue Book.

Sections:
- P-31: Modification of Degree of Hydrogenation
- P-32: Prefixes for Substituent Groups
- P-33: Suffixes
- P-34: Functional Parent Compounds
- P-35: Prefixes Corresponding to Characteristic Groups
"""

from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


# =============================================================================
# P-31 MODIFICATION OF DEGREE OF HYDROGENATION
# =============================================================================

def format_unsaturation(double_bonds: List[int], triple_bonds: List[int]) -> Tuple[str, str]:
    """
    Format unsaturation suffixes (P-31.1).

    Returns:
        Tuple of (locants_string, suffix_string)

    Examples:
        ([2], []) -> ("2", "ene")
        ([2, 4], []) -> ("2,4", "diene")
        ([2], [4]) -> ("2", "en"), ("4", "yne")
    """
    from .p1_general import get_multiplicative_prefix, format_locants

    parts = []

    # Double bonds
    if double_bonds:
        mult = get_multiplicative_prefix(len(double_bonds))
        locants = format_locants(double_bonds)
        suffix = f"{mult}en" if len(double_bonds) > 1 else "en"
        parts.append((locants, suffix))

    # Triple bonds
    if triple_bonds:
        mult = get_multiplicative_prefix(len(triple_bonds))
        locants = format_locants(triple_bonds)
        suffix = f"{mult}yn" if len(triple_bonds) > 1 else "yn"
        parts.append((locants, suffix))

    return parts


# P-31.2 Hydro/dehydro prefixes
def format_hydro_prefix(positions: List[int], n_hydrogens: int = 2) -> str:
    """
    Format 'hydro' prefix for saturation (P-31.2.3).

    Example: 1,2-dihydro, 1,2,3,4-tetrahydro
    """
    from .p1_general import get_multiplicative_prefix, format_locants

    mult = get_multiplicative_prefix(len(positions) * (n_hydrogens // 2))
    locants = format_locants(positions)

    return f"{locants}-{mult}hydro"


def format_dehydro_prefix(positions: List[int]) -> str:
    """
    Format 'dehydro' prefix for introduction of unsaturation (P-31.2.4).
    """
    from .p1_general import get_multiplicative_prefix, format_locants

    mult = get_multiplicative_prefix(len(positions))
    locants = format_locants(positions)

    return f"{locants}-{mult}dehydro"


# =============================================================================
# P-33 SUFFIXES
# =============================================================================

# P-33.2.1 Basic functional suffixes (in order of seniority)
FUNCTIONAL_SUFFIXES = {
    # Acids (highest priority in this group)
    'carboxylic_acid': {'suffix': 'oic acid', 'prefix': 'carboxy', 'detachable': True},
    'sulfonic_acid': {'suffix': 'sulfonic acid', 'prefix': 'sulfo', 'detachable': True},
    'sulfinic_acid': {'suffix': 'sulfinic acid', 'prefix': 'sulfino', 'detachable': True},

    # Acid derivatives
    'acid_anhydride': {'suffix': 'oic anhydride', 'detachable': True},
    'ester': {'suffix': 'oate', 'detachable': True},
    'acid_halide': {'suffix': 'oyl halide', 'detachable': True},
    'amide': {'suffix': 'amide', 'prefix': 'carbamoyl', 'detachable': True},
    'hydrazide': {'suffix': 'ohydrazide', 'detachable': True},
    'nitrile': {'suffix': 'nitrile', 'prefix': 'cyano', 'detachable': False},

    # Aldehydes and ketones
    'aldehyde': {'suffix': 'al', 'prefix': 'formyl/oxo', 'detachable': False},
    'ketone': {'suffix': 'one', 'prefix': 'oxo', 'detachable': False},

    # Alcohols and thiols
    'alcohol': {'suffix': 'ol', 'prefix': 'hydroxy', 'detachable': False},
    'thiol': {'suffix': 'thiol', 'prefix': 'sulfanyl/mercapto', 'detachable': False},

    # Amines
    'amine': {'suffix': 'amine', 'prefix': 'amino', 'detachable': False},
    'imine': {'suffix': 'imine', 'prefix': 'imino', 'detachable': False},
}


# Suffix modifications for carbon chains vs rings
def get_suffix_for_chain(functional_group: str, chain_length: int) -> str:
    """
    Get appropriate suffix for a carbon chain with functional group.

    Example: methanol, ethanoic acid, propanone
    """
    data = FUNCTIONAL_SUFFIXES.get(functional_group, {})
    suffix = data.get('suffix', '')

    if functional_group == 'carboxylic_acid':
        return 'oic acid'
    elif functional_group == 'aldehyde':
        return 'al'
    elif functional_group == 'ketone':
        return 'one'
    elif functional_group == 'alcohol':
        return 'ol'
    elif functional_group == 'amine':
        return 'amine'
    elif functional_group == 'nitrile':
        return 'nitrile'

    return suffix


# =============================================================================
# P-34 FUNCTIONAL PARENT COMPOUNDS
# =============================================================================

# P-34.1 Retained functional parent compounds
RETAINED_FUNCTIONAL_PARENTS = {
    # Carboxylic acids
    'formic acid': {'formula': 'HCOOH', 'systematic': 'methanoic acid'},
    'acetic acid': {'formula': 'CH3COOH', 'systematic': 'ethanoic acid'},
    'propionic acid': {'formula': 'C2H5COOH', 'systematic': 'propanoic acid'},
    'butyric acid': {'formula': 'C3H7COOH', 'systematic': 'butanoic acid'},
    'valeric acid': {'formula': 'C4H9COOH', 'systematic': 'pentanoic acid'},
    'caproic acid': {'formula': 'C5H11COOH', 'systematic': 'hexanoic acid'},
    'benzoic acid': {'formula': 'C6H5COOH'},
    'oxalic acid': {'formula': 'HOOC-COOH', 'systematic': 'ethanedioic acid'},
    'malonic acid': {'formula': 'HOOC-CH2-COOH', 'systematic': 'propanedioic acid'},
    'succinic acid': {'formula': 'HOOC-(CH2)2-COOH', 'systematic': 'butanedioic acid'},
    'glutaric acid': {'formula': 'HOOC-(CH2)3-COOH', 'systematic': 'pentanedioic acid'},
    'adipic acid': {'formula': 'HOOC-(CH2)4-COOH', 'systematic': 'hexanedioic acid'},

    # Aldehydes
    'formaldehyde': {'formula': 'HCHO', 'systematic': 'methanal'},
    'acetaldehyde': {'formula': 'CH3CHO', 'systematic': 'ethanal'},
    'benzaldehyde': {'formula': 'C6H5CHO'},

    # Ketones
    'acetone': {'formula': '(CH3)2CO', 'systematic': 'propan-2-one'},
    'acetophenone': {'formula': 'C6H5COCH3', 'systematic': '1-phenylethan-1-one'},
    'benzophenone': {'formula': '(C6H5)2CO', 'systematic': 'diphenylmethanone'},

    # Alcohols
    'methanol': {'formula': 'CH3OH'},
    'ethanol': {'formula': 'C2H5OH'},
    'phenol': {'formula': 'C6H5OH'},
    'glycerol': {'formula': 'HOCH2CH(OH)CH2OH', 'systematic': 'propane-1,2,3-triol'},

    # Amines
    'aniline': {'formula': 'C6H5NH2', 'systematic': 'phenylamine'},

    # Others
    'urea': {'formula': '(NH2)2CO', 'systematic': 'carbonyl diamide'},
    'acetamide': {'formula': 'CH3CONH2', 'systematic': 'ethanamide'},
}


# =============================================================================
# P-35 PREFIXES FOR CHARACTERISTIC GROUPS
# =============================================================================

# P-35.2.1 Retained traditional prefixes
CHARACTERISTIC_GROUP_PREFIXES = {
    # Oxygen groups
    'hydroxy': {'group': '-OH', 'suffix': 'ol'},
    'oxo': {'group': '=O', 'suffix': 'one/al'},
    'carboxy': {'group': '-COOH', 'suffix': 'oic acid'},
    'formyl': {'group': '-CHO', 'suffix': 'al'},
    'acetyl': {'group': '-COCH3'},
    'benzoyl': {'group': '-COC6H5'},

    # Nitrogen groups
    'amino': {'group': '-NH2', 'suffix': 'amine'},
    'imino': {'group': '=NH', 'suffix': 'imine'},
    'nitro': {'group': '-NO2'},
    'nitroso': {'group': '-NO'},
    'cyano': {'group': '-CN', 'suffix': 'nitrile'},
    'isocyano': {'group': '-NC'},
    'azido': {'group': '-N3'},
    'diazo': {'group': '=N2'},

    # Sulfur groups
    'sulfanyl': {'group': '-SH', 'suffix': 'thiol', 'alt': 'mercapto'},
    'sulfinyl': {'group': '-S(O)-'},
    'sulfonyl': {'group': '-S(O)2-'},
    'sulfo': {'group': '-SO3H', 'suffix': 'sulfonic acid'},

    # Halogen groups
    'fluoro': {'group': '-F'},
    'chloro': {'group': '-Cl'},
    'bromo': {'group': '-Br'},
    'iodo': {'group': '-I'},

    # Other groups
    'phosphono': {'group': '-PO(OH)2'},
    'silyl': {'group': '-SiH3'},
}


# P-35.2.2 Substituents from parent hydrides
HYDRIDE_DERIVED_PREFIXES = {
    # Alkyl groups
    'methyl': {'from': 'methane', 'formula': '-CH3'},
    'ethyl': {'from': 'ethane', 'formula': '-C2H5'},
    'propyl': {'from': 'propane', 'formula': '-C3H7'},
    'butyl': {'from': 'butane', 'formula': '-C4H9'},
    'pentyl': {'from': 'pentane', 'formula': '-C5H11'},
    'hexyl': {'from': 'hexane', 'formula': '-C6H13'},

    # Branched alkyl (retained names)
    'isopropyl': {'from': 'propane', 'systematic': 'propan-2-yl'},
    'isobutyl': {'from': 'butane', 'systematic': '2-methylpropyl'},
    'sec-butyl': {'from': 'butane', 'systematic': 'butan-2-yl'},
    'tert-butyl': {'from': 'butane', 'systematic': '2-methylpropan-2-yl'},
    'neopentyl': {'from': 'pentane', 'systematic': '2,2-dimethylpropyl'},

    # Unsaturated
    'vinyl': {'from': 'ethene', 'systematic': 'ethenyl'},
    'allyl': {'from': 'propene', 'systematic': 'prop-2-en-1-yl'},
    'ethynyl': {'from': 'ethyne'},
    'propargyl': {'from': 'propyne', 'systematic': 'prop-2-yn-1-yl'},

    # Aryl
    'phenyl': {'from': 'benzene', 'formula': '-C6H5'},
    'benzyl': {'from': 'toluene', 'systematic': 'phenylmethyl'},
    'naphthyl': {'from': 'naphthalene'},
    'tolyl': {'from': 'toluene', 'systematic': 'methylphenyl'},
}


def get_substituent_prefix(name: str) -> Optional[str]:
    """
    Get the prefix form of a substituent.

    Args:
        name: Name of the substituent

    Returns:
        Prefix string or None if not found
    """
    if name in CHARACTERISTIC_GROUP_PREFIXES:
        return name
    if name in HYDRIDE_DERIVED_PREFIXES:
        return name
    return None


def get_suffix_for_group(group_type: str) -> Optional[str]:
    """
    Get the suffix for a characteristic group type.
    """
    if group_type in FUNCTIONAL_SUFFIXES:
        return FUNCTIONAL_SUFFIXES[group_type].get('suffix')
    return None


# =============================================================================
# COMPOUND AND COMPLEX SUBSTITUENTS (P-35.3, P-35.4)
# =============================================================================

def build_compound_substituent(base: str, substituents: List[Tuple[int, str]]) -> str:
    """
    Build a compound substituent name (P-35.3).

    Args:
        base: Base substituent name
        substituents: List of (position, substituent_name) tuples

    Example:
        ('methyl', [(2, 'chloro')]) -> '2-chloromethyl'
    """
    from .p1_general import format_locants

    # Sort substituents alphabetically
    sorted_subs = sorted(substituents, key=lambda x: x[1])

    parts = []
    for pos, sub in sorted_subs:
        parts.append(f"{pos}-{sub}")

    return ''.join(parts) + base