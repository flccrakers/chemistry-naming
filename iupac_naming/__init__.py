"""
IUPAC Nomenclature for Organic Compounds
=========================================

Implementation of IUPAC 2013 Blue Book nomenclature rules.

Structure follows the Blue Book chapters:
- P-1: General Principles, Rules, and Conventions
- P-2: Parent Hydrides
- P-3: Characteristic (Functional) and Substituent Groups
- P-4: Rules for Name Construction
- P-5: Selecting Preferred IUPAC Names
- P-6: Applications to Specific Classes of Compounds
- P-7: Radicals, Ions, and Related Species
- P-8: Isotopically Modified Compounds
- P-9: Specification of Configuration and Conformation

Usage:
    from iupac_nomenclature import IUPACNamer
    
    namer = IUPACNamer()
    name = namer.name_from_smiles('CCO')  # -> ethanol

Author: Generated with Claude AI
Reference: IUPAC Recommendations 2013 (Blue Book)
"""

__version__ = "0.1.0"
__author__ = "IUPAC Nomenclature Project"

from .namer import IUPACNamer
from .data_structures import (
    AtomInfo, RingInfo, FusedRingSystem, SpiroSystem,
    CharacteristicGroup, MoleculeAnalysis
)

__all__ = [
    "IUPACNamer",
    "AtomInfo",
    "RingInfo", 
    "FusedRingSystem",
    "SpiroSystem",
    "CharacteristicGroup",
    "MoleculeAnalysis",
]
