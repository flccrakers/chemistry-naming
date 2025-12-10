"""
IUPAC Nomenclature Rules
========================

Implementation of IUPAC 2013 Blue Book chapters P-1 through P-9.
"""

from .p1_general import (
    NameType, NomenclatureOperation, NomenclatureType,
    get_standard_bonding_number, format_lambda_convention,
    get_multiplicative_prefix, format_locants, compare_locant_sets,
    alphanumerical_sort_key, format_indicated_hydrogen,
    NameWriter, add_primes
)

from .p2_parent_hydrides import (
    MONONUCLEAR_HYDRIDES, RETAINED_HETEROCYCLES,
    RETAINED_FUSED_HYDROCARBONS, RETAINED_FUSED_HETEROCYCLES,
    FUSION_PREFIXES, BIVALENT_BRIDGES, RING_ASSEMBLY_PREFIXES,
    HANTZSCH_WIDMAN_PREFIXES, HW_STEMS_SATURATED, HW_STEMS_UNSATURATED,
    get_acyclic_hydride_name, get_cycloalkane_name,
    get_hantzsch_widman_name, get_ring_assembly_name
)

from .p3_characteristic_groups import (
    FUNCTIONAL_SUFFIXES, RETAINED_FUNCTIONAL_PARENTS,
    CHARACTERISTIC_GROUP_PREFIXES, HYDRIDE_DERIVED_PREFIXES,
    format_unsaturation, format_hydro_prefix, format_dehydro_prefix,
    get_substituent_prefix, get_suffix_for_group,
    build_compound_substituent
)

from .p4_name_construction import (
    CLASS_SENIORITY, get_class_seniority, compare_class_seniority,
    compare_ring_seniority, compare_fused_system_seniority,
    select_principal_chain, apply_low_locant_rule,
    select_preferred_name, select_principal_substituent_chain,
    determine_parent_structure
)

__all__ = [
    # P-1 General
    'NameType', 'NomenclatureOperation', 'NomenclatureType',
    'get_standard_bonding_number', 'format_lambda_convention',
    'get_multiplicative_prefix', 'format_locants', 'compare_locant_sets',
    'alphanumerical_sort_key', 'format_indicated_hydrogen',
    'NameWriter', 'add_primes',
    
    # P-2 Parent Hydrides
    'MONONUCLEAR_HYDRIDES', 'RETAINED_HETEROCYCLES',
    'RETAINED_FUSED_HYDROCARBONS', 'RETAINED_FUSED_HETEROCYCLES',
    'FUSION_PREFIXES', 'BIVALENT_BRIDGES', 'RING_ASSEMBLY_PREFIXES',
    'HANTZSCH_WIDMAN_PREFIXES', 'HW_STEMS_SATURATED', 'HW_STEMS_UNSATURATED',
    'get_acyclic_hydride_name', 'get_cycloalkane_name',
    'get_hantzsch_widman_name', 'get_ring_assembly_name',
    
    # P-3 Characteristic Groups
    'FUNCTIONAL_SUFFIXES', 'RETAINED_FUNCTIONAL_PARENTS',
    'CHARACTERISTIC_GROUP_PREFIXES', 'HYDRIDE_DERIVED_PREFIXES',
    'format_unsaturation', 'format_hydro_prefix', 'format_dehydro_prefix',
    'get_substituent_prefix', 'get_suffix_for_group',
    'build_compound_substituent',
    
    # P-4 Name Construction
    'CLASS_SENIORITY', 'get_class_seniority', 'compare_class_seniority',
    'compare_ring_seniority', 'compare_fused_system_seniority',
    'select_principal_chain', 'apply_low_locant_rule',
    'select_preferred_name', 'select_principal_substituent_chain',
    'determine_parent_structure',
]
