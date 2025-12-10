"""
IUPAC Nomenclature Utilities
============================

Helper functions for IUPAC nomenclature generation.
"""

from typing import List, Dict, Set, Tuple, Optional
import re


def normalize_locants(locants: List[str]) -> List[str]:
    """
    Normalize locant format for comparison.
    
    Handles:
    - Numeric locants (1, 2, 3...)
    - Primed locants (1', 2'', 3'''...)
    - Letter locants (a, b, c...)
    
    Args:
        locants: List of locant strings
    
    Returns:
        Sorted, normalized locants
    """
    def locant_key(loc: str) -> tuple:
        # Extract numeric part and prime count
        match = re.match(r"(\d+)('+)?([a-z])?", loc)
        if match:
            num = int(match.group(1))
            primes = len(match.group(2)) if match.group(2) else 0
            letter = match.group(3) or ''
            return (primes, num, letter)
        
        # Letter-only locant
        if loc.isalpha():
            return (0, ord(loc), '')
        
        return (999, 0, loc)
    
    return sorted(locants, key=locant_key)


def compare_locant_sequences(seq1: List[str], seq2: List[str]) -> int:
    """
    Compare two sequences of locants for the low locant rule.
    
    Reference: P-14.3.5
    
    Returns:
        -1 if seq1 < seq2 (seq1 is lower)
         0 if equal
         1 if seq1 > seq2 (seq2 is lower)
    """
    # Normalize both sequences
    norm1 = normalize_locants(seq1)
    norm2 = normalize_locants(seq2)
    
    # Compare position by position
    for l1, l2 in zip(norm1, norm2):
        k1 = _locant_sort_key(l1)
        k2 = _locant_sort_key(l2)
        
        if k1 < k2:
            return -1
        elif k1 > k2:
            return 1
    
    # If all compared locants equal, shorter sequence is lower
    if len(norm1) < len(norm2):
        return -1
    elif len(norm1) > len(norm2):
        return 1
    
    return 0


def _locant_sort_key(loc: str) -> tuple:
    """Get sort key for a single locant."""
    match = re.match(r"(\d+)('+)?", loc)
    if match:
        num = int(match.group(1))
        primes = len(match.group(2)) if match.group(2) else 0
        return (primes, num)
    return (999, ord(loc[0]) if loc else 999)


def format_locant_list(locants: List[str], separator: str = ',') -> str:
    """
    Format a list of locants as a string.
    
    Args:
        locants: List of locant strings
        separator: Separator character (default comma)
    
    Returns:
        Formatted string like "1,2,3" or "1',2',3'"
    """
    if not locants:
        return ""
    
    normalized = normalize_locants(locants)
    return separator.join(normalized)


def elide_vowel(word: str, suffix: str) -> str:
    """
    Apply vowel elision rules for IUPAC naming.
    
    Reference: P-16.7
    
    Rules:
    - Final 'a' of a name is elided before vowel
    - Final 'e' of parent name elided before suffix starting with 'a', 'e', 'i', 'o', 'y'
    - 'o' of 'methano', 'ethano' etc. elided before 'a' or 'o'
    
    Args:
        word: Base word
        suffix: Suffix to add
    
    Returns:
        Combined word with appropriate elision
    """
    if not word or not suffix:
        return word + suffix
    
    vowels = 'aeiou'
    suffix_vowels = 'aeioy'
    
    # Rule: final 'a' elided before vowel
    if word.endswith('a') and suffix[0] in vowels:
        return word[:-1] + suffix
    
    # Rule: final 'e' elided before certain vowels
    if word.endswith('e') and suffix[0] in suffix_vowels:
        return word[:-1] + suffix
    
    # Rule: 'o' in bridge prefixes elided before 'a' or 'o'
    bridge_prefixes = ['methano', 'ethano', 'propano', 'butano', 'epoxy', 'epithio']
    for prefix in bridge_prefixes:
        if word.endswith(prefix) and suffix[0] in 'ao':
            return word[:-1] + suffix
    
    return word + suffix


def split_name_components(name: str) -> Dict[str, str]:
    """
    Split an IUPAC name into its components.
    
    Useful for parsing and analyzing existing names.
    
    Args:
        name: IUPAC name string
    
    Returns:
        Dictionary with keys: prefixes, parent, suffixes, locants
    """
    components = {
        'prefixes': [],
        'parent': '',
        'suffixes': [],
        'locants': [],
        'original': name
    }
    
    # Extract locants (numbers and primes at start or before hyphens)
    locant_pattern = r"(\d+'*)"
    locants = re.findall(locant_pattern, name)
    components['locants'] = locants
    
    # Remove locants for further parsing
    cleaned = re.sub(r"\d+'*[,-]?", '', name)
    
    # Common parent hydrides
    parents = [
        'benzene', 'naphthalene', 'anthracene', 'phenanthrene',
        'pyridine', 'pyrimidine', 'pyrazine', 'pyridazine',
        'furan', 'thiophene', 'pyrrole',
        'imidazole', 'pyrazole', 'oxazole', 'thiazole',
        'quinoline', 'isoquinoline', 'indole',
        'methane', 'ethane', 'propane', 'butane', 'pentane',
        'hexane', 'heptane', 'octane', 'nonane', 'decane',
        'methanol', 'ethanol',
        'phenazine', 'acridine', 'carbazole',
    ]
    
    for parent in sorted(parents, key=len, reverse=True):
        if parent in cleaned.lower():
            components['parent'] = parent
            break
    
    return components


def get_element_valence(symbol: str) -> int:
    """
    Get standard valence for an element.
    
    Reference: P-14.1
    
    Args:
        symbol: Element symbol (e.g., 'C', 'N', 'O')
    
    Returns:
        Standard valence
    """
    valences = {
        'H': 1,
        'B': 3,
        'C': 4,
        'N': 3,
        'O': 2,
        'F': 1,
        'Si': 4,
        'P': 3,
        'S': 2,
        'Cl': 1,
        'Br': 1,
        'I': 1,
        'Se': 2,
        'Te': 2,
        'As': 3,
        'Sb': 3,
        'Bi': 3,
    }
    return valences.get(symbol, 4)


def get_heteroatom_priority(symbol: str) -> int:
    """
    Get priority value for a heteroatom.
    
    Reference: P-44.2.2 (seniority order)
    
    Higher value = higher priority (appears first in name)
    
    Order: O > S > Se > Te > N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B
    
    Args:
        symbol: Element symbol
    
    Returns:
        Priority value (higher = more senior)
    """
    priorities = {
        'O': 100,
        'S': 90,
        'Se': 80,
        'Te': 70,
        'N': 60,
        'P': 50,
        'As': 40,
        'Sb': 35,
        'Bi': 30,
        'Si': 25,
        'Ge': 20,
        'Sn': 15,
        'Pb': 10,
        'B': 5,
        'C': 0,  # Carbon has lowest priority (not really a heteroatom)
    }
    return priorities.get(symbol, 0)


def is_aromatic_ring(elements: tuple, bond_types: tuple) -> bool:
    """
    Determine if a ring is aromatic based on elements and bond types.
    
    Simple heuristic - real aromaticity is complex.
    
    Args:
        elements: Tuple of element symbols in ring
        bond_types: Tuple of bond type integers
    
    Returns:
        True if likely aromatic
    """
    # Must have alternating bonds or be benzene-like
    size = len(elements)
    
    # 6-membered all carbons with alternating bonds
    if size == 6:
        if all(e == 'C' for e in elements):
            return True
        # May have N, O, S
        if sum(1 for e in elements if e in 'NOS') <= 2:
            return True
    
    # 5-membered with heteroatom
    if size == 5:
        hetero_count = sum(1 for e in elements if e in 'NOS')
        if hetero_count >= 1:
            return True
    
    return False


def generate_primed_locants(base: int, prime_level: int) -> str:
    """
    Generate locant with primes.
    
    Args:
        base: Numeric locant
        prime_level: Number of primes (0 = none, 1 = ', 2 = '', etc.)
    
    Returns:
        Locant string like "1", "1'", "1''", etc.
    """
    if prime_level == 0:
        return str(base)
    return str(base) + "'" * prime_level


def parse_fusion_locants(locant_str: str) -> Dict:
    """
    Parse fusion locant notation.
    
    Examples:
    - [2,3-b] -> {'attached': [2,3], 'base': 'b'}
    - [4',5':5,6] -> {'attached': [4',5'], 'base': [5,6]}
    
    Args:
        locant_str: Fusion locant string
    
    Returns:
        Dictionary with parsed components
    """
    result = {
        'attached_locants': [],
        'base_locants': [],
        'face_letter': None,
        'raw': locant_str
    }
    
    # Remove brackets
    inner = locant_str.strip('[]')
    
    # Check for colon (full locants) or hyphen-letter (face notation)
    if ':' in inner:
        # Full notation: 4',5':5,6
        parts = inner.split(':')
        result['attached_locants'] = parts[0].split(',')
        result['base_locants'] = parts[1].split(',')
    elif '-' in inner:
        # Face notation: 2,3-b
        parts = inner.split('-')
        result['attached_locants'] = parts[0].split(',')
        result['face_letter'] = parts[1] if len(parts) > 1 else None
    else:
        # Just locants
        result['attached_locants'] = inner.split(',')
    
    return result
