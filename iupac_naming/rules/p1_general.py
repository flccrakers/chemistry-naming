"""
P-1 GENERAL PRINCIPLES, RULES, AND CONVENTIONS
==============================================

This module implements the fundamental rules from Chapter P-1 of the
IUPAC 2013 Blue Book.

Sections:
- P-10: Introduction
- P-11: Scope of nomenclature for organic compounds
- P-12: Preferred, Preselected, and Retained IUPAC Names
- P-13: Operations in Nomenclature of Organic Compounds
- P-14: General Rules
- P-15: Types of Nomenclature
- P-16: Name Writing
"""

from typing import List, Dict, Optional, Union, Tuple
from dataclasses import dataclass
from enum import Enum, auto


# =============================================================================
# P-12 PREFERRED, PRESELECTED, AND RETAINED IUPAC NAMES
# =============================================================================

class NameType(Enum):
    """
    Types of IUPAC names (P-12).
    
    P-12.1: Preferred IUPAC Names (PINs) - unique names for regulatory purposes
    P-12.2: Preselected Names - names for parent structures  
    P-12.3: Retained Names - traditional names kept for specific compounds
    """
    PREFERRED = auto()      # PIN
    PRESELECTED = auto()    # For parent structures
    RETAINED = auto()       # Traditional names
    ACCEPTABLE = auto()     # Valid but not preferred
    GENERAL = auto()        # For general nomenclature


# =============================================================================
# P-13 NOMENCLATURE OPERATIONS
# =============================================================================

class NomenclatureOperation(Enum):
    """
    Operations used in systematic nomenclature (P-13).
    """
    # P-13.1 Substitutive operation
    SUBSTITUTIVE = auto()       # Replacing H atoms
    
    # P-13.2 Replacement operation (skeletal)
    REPLACEMENT = auto()        # 'a' nomenclature (oxa, aza, etc.)
    
    # P-13.3 Additive operation
    ADDITIVE_PREFIX = auto()    # e.g., epoxy, peroxy
    ADDITIVE_SUFFIX = auto()    # e.g., -oxide
    ADDITIVE_WORD = auto()      # Separate word additions
    
    # P-13.4 Subtractive operation
    SUBTRACTIVE_SUFFIX = auto() # e.g., -ene, -yne
    SUBTRACTIVE_CHANGE = auto() # Change in ending
    SUBTRACTIVE_PREFIX = auto() # dehydro, nor
    
    # P-13.5 Conjunctive operation
    CONJUNCTIVE = auto()        # Joining ring and chain
    
    # P-13.6 Multiplicative operation
    MULTIPLICATIVE = auto()     # Identical units linked by central atom
    
    # P-13.7 Fusion operation
    FUSION = auto()             # Fused ring names


# =============================================================================
# P-14.1 BONDING NUMBERS
# =============================================================================

# Standard bonding numbers (P-14.1.2)
STANDARD_BONDING_NUMBERS = {
    'H': 1, 'B': 3, 'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
    'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,
    'O': 2, 'S': 2, 'Se': 2, 'Te': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1
}


def get_standard_bonding_number(element: str) -> int:
    """
    Get the standard bonding number for an element (P-14.1.2).
    """
    return STANDARD_BONDING_NUMBERS.get(element, 4)


def format_lambda_convention(element: str, bonding_number: int) -> str:
    """
    Format nonstandard bonding number using λ-convention (P-14.1.3).
    
    Example: λ5-phosphane for P with bonding number 5
    """
    standard = get_standard_bonding_number(element)
    if bonding_number != standard:
        return f"λ{bonding_number}-"
    return ""


# =============================================================================
# P-14.2 MULTIPLICATIVE PREFIXES
# =============================================================================

def get_multiplicative_prefix(n: int, use_bis: bool = False) -> str:
    """
    Get multiplicative prefix for a number (P-14.2).
    
    Args:
        n: Number of items
        use_bis: If True, use bis/tris/tetrakis form (P-14.2.3)
    
    Examples:
        2 -> "di" or "bis"
        3 -> "tri" or "tris"
        4 -> "tetra" or "tetrakis"
    """
    if n == 1:
        return ""
    
    if use_bis:
        bis_prefixes = {
            2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis",
            6: "hexakis", 7: "heptakis", 8: "octakis", 9: "nonakis",
            10: "decakis", 11: "undecakis", 12: "dodecakis"
        }
        return bis_prefixes.get(n, f"{n}kis")
    
    basic_prefixes = {
        2: "di", 3: "tri", 4: "tetra", 5: "penta",
        6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca",
        11: "undeca", 12: "dodeca", 13: "trideca", 14: "tetradeca",
        15: "pentadeca", 16: "hexadeca", 17: "heptadeca", 18: "octadeca",
        19: "nonadeca", 20: "icosa", 21: "henicosa", 22: "docosa",
        23: "tricosa", 30: "triaconta", 40: "tetraconta", 50: "pentaconta",
        60: "hexaconta", 70: "heptaconta", 80: "octaconta", 90: "nonaconta",
        100: "hecta"
    }
    
    if n <= 22:
        return basic_prefixes.get(n, f"{n}")
    elif n < 100:
        tens = (n // 10) * 10
        units = n % 10
        tens_prefix = basic_prefixes.get(tens, "")
        units_prefix = basic_prefixes.get(units, "")
        if units == 1:
            return f"hen{tens_prefix}"
        elif units == 2:
            return f"do{tens_prefix}"
        else:
            return f"{units_prefix}{tens_prefix}"
    
    return str(n)


# =============================================================================
# P-14.3 LOCANTS
# =============================================================================

def format_locants(locants: List[Union[int, str]], separator: str = ",") -> str:
    """
    Format a list of locants (P-14.3).
    
    Rules:
    - P-14.3.2: Locants immediately precede the part of the name to which they relate
    - P-14.3.3: Citation of locants
    
    Examples:
        [1, 2, 3] -> "1,2,3"
        [1, "2a", 3] -> "1,2a,3"
    """
    return separator.join(str(l) for l in locants)


def compare_locant_sets(set1: List[Union[int, str]], 
                        set2: List[Union[int, str]]) -> int:
    """
    Compare two sets of locants to determine which is lower (P-14.3.5).
    
    Returns:
        -1 if set1 is lower (preferred)
        0 if equal
        1 if set2 is lower (preferred)
    """
    # Convert to comparable format
    def locant_key(loc):
        if isinstance(loc, int):
            return (loc, "")
        # Handle primed locants like 1', 2''
        base = ""
        primes = 0
        letters = ""
        for c in str(loc):
            if c.isdigit():
                base += c
            elif c == "'":
                primes += 1
            else:
                letters += c
        return (int(base) if base else 0, primes, letters)
    
    sorted1 = sorted(set1, key=locant_key)
    sorted2 = sorted(set2, key=locant_key)
    
    for l1, l2 in zip(sorted1, sorted2):
        k1, k2 = locant_key(l1), locant_key(l2)
        if k1 < k2:
            return -1
        elif k1 > k2:
            return 1
    
    # If all compared are equal, shorter set wins
    if len(set1) < len(set2):
        return -1
    elif len(set1) > len(set2):
        return 1
    
    return 0


# =============================================================================
# P-14.5 ALPHANUMERICAL ORDER
# =============================================================================

def alphanumerical_sort_key(name: str) -> str:
    """
    Generate a sort key for alphanumerical ordering (P-14.5).
    
    Rules:
    - P-14.5.1: Simple prefixes are alphabetized based on the complete prefix
    - P-14.5.2: Multiplicative prefixes (di, tri, etc.) are NOT considered
    - P-14.5.3: Locants, Greek letters, etc. are ignored
    """
    # Remove multiplicative prefixes for sorting
    prefixes_to_remove = [
        'di', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona', 'deca',
        'bis', 'tris', 'tetrakis', 'pentakis', 'hexakis'
    ]
    
    result = name.lower()
    
    # Remove locants (digits and commas at start)
    import re
    result = re.sub(r'^[\d,\-\']+', '', result)
    
    # Remove multiplicative prefixes
    for prefix in sorted(prefixes_to_remove, key=len, reverse=True):
        if result.startswith(prefix):
            result = result[len(prefix):]
            break
    
    return result


# =============================================================================
# P-14.7 INDICATED HYDROGEN
# =============================================================================

def format_indicated_hydrogen(position: Union[int, str]) -> str:
    """
    Format indicated hydrogen notation (P-14.7.1).
    
    Example: 1H, 2H, 4aH
    """
    return f"{position}H"


def format_added_indicated_hydrogen(position: Union[int, str]) -> str:
    """
    Format added indicated hydrogen notation (P-14.7.2).
    
    Example: 9aH (in brackets when needed)
    """
    return f"{position}H"


# =============================================================================
# P-15 TYPES OF NOMENCLATURE
# =============================================================================

class NomenclatureType(Enum):
    """
    Types of IUPAC nomenclature (P-15).
    """
    # P-15.1 Substitutive nomenclature (most common)
    SUBSTITUTIVE = auto()
    
    # P-15.2 Functional class nomenclature
    FUNCTIONAL_CLASS = auto()
    
    # P-15.3 Multiplicative nomenclature
    MULTIPLICATIVE = auto()
    
    # P-15.4 Skeletal replacement ('a') nomenclature
    REPLACEMENT = auto()
    
    # P-15.5 Functional replacement nomenclature
    FUNCTIONAL_REPLACEMENT = auto()
    
    # P-15.6 Conjunctive nomenclature
    CONJUNCTIVE = auto()


# =============================================================================
# P-16 NAME WRITING
# =============================================================================

class NameWriter:
    """
    Utilities for writing IUPAC names correctly (P-16).
    """
    
    @staticmethod
    def format_locant_with_hyphen(locant: Union[int, str], term: str) -> str:
        """
        P-16.2.4: Format locant-term with proper hyphen.
        """
        return f"{locant}-{term}"
    
    @staticmethod
    def join_with_comma(items: List[str]) -> str:
        """
        P-16.2.1: Join items with commas.
        """
        return ",".join(items)
    
    @staticmethod
    def needs_parentheses(name: str) -> bool:
        """
        P-16.5.1: Determine if parentheses are needed.
        
        Parentheses are used for:
        - Compound/complex substituents
        - When needed to avoid ambiguity
        """
        # Complex substituents containing locants or hyphens
        if '-' in name and any(c.isdigit() for c in name):
            return True
        return False
    
    @staticmethod
    def elide_vowel(prefix: str, suffix: str) -> str:
        """
        P-16.7: Vowel elision rules.
        
        Terminal 'a' is elided before 'a', 'e', 'i', 'o', 'u', 'y'
        Terminal 'o' is elided before 'a', 'o'
        """
        if not prefix or not suffix:
            return prefix + suffix
        
        vowels_for_a = set('aeiouy')
        vowels_for_o = set('ao')
        
        if prefix[-1] == 'a' and suffix[0].lower() in vowels_for_a:
            return prefix[:-1] + suffix
        elif prefix[-1] == 'o' and suffix[0].lower() in vowels_for_o:
            return prefix[:-1] + suffix
        
        return prefix + suffix


# =============================================================================
# P-16.9 PRIMES
# =============================================================================

def add_primes(locant: Union[int, str], level: int = 1) -> str:
    """
    Add primes to a locant (P-16.9).
    
    Args:
        locant: The base locant
        level: Number of primes to add
    
    Examples:
        add_primes(1, 1) -> "1'"
        add_primes(1, 2) -> "1''"
        add_primes("2a", 1) -> "2a'"
    """
    primes = "'" * level
    return f"{locant}{primes}"
