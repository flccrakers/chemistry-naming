"""
IUPAC Nomenclature Analyzers
============================

Modules for analyzing molecular structures.
"""

from .molecular_analyzer import MolecularAnalyzer
from .component_matcher import ComponentMatcher
from .peripheral_numbering import (
    PeripheralNumberer,
    compute_fusion_descriptor,
    format_fusion_name,
    order_components_by_distance
)

__all__ = [
    'MolecularAnalyzer',
    'ComponentMatcher',
    'PeripheralNumberer',
    'compute_fusion_descriptor',
    'format_fusion_name',
    'order_components_by_distance',
]