"""
Common data structures for IUPAC nomenclature.

These dataclasses represent molecular structures and their components
as needed for nomenclature generation.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set, Union, FrozenSet
from enum import Enum, auto


# =============================================================================
# P-14.1 BONDING NUMBERS
# =============================================================================

class BondingNumber(Enum):
    """
    Standard and nonstandard bonding numbers.
    Reference: P-14.1
    """
    # Standard bonding numbers (P-14.1.2)
    STANDARD = auto()
    # Nonstandard bonding numbers using λ-convention (P-14.1.3)
    LAMBDA = auto()


# =============================================================================
# P-14.2 MULTIPLICATIVE PREFIXES
# =============================================================================

MULTIPLICATIVE_PREFIXES = {
    # P-14.2.1 Basic multiplicative prefixes
    1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",
    6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca",
    11: "undeca", 12: "dodeca", 13: "trideca", 14: "tetradeca",
    15: "pentadeca", 16: "hexadeca", 17: "heptadeca", 18: "octadeca",
    19: "nonadeca", 20: "icosa", 21: "henicosa", 22: "docosa",
    30: "triaconta", 40: "tetraconta", 50: "pentaconta",
    100: "hecta", 200: "dicta", 1000: "kilia"
}

MULTIPLICATIVE_PREFIXES_BIS = {
    # P-14.2.1 For complex groups
    2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis",
    6: "hexakis", 7: "heptakis", 8: "octakis", 9: "nonakis",
    10: "decakis"
}


# =============================================================================
# P-21 MONONUCLEAR AND ACYCLIC PARENT HYDRIDES
# =============================================================================

@dataclass
class AtomInfo:
    """Information about a single atom in the molecule."""
    idx: int
    symbol: str
    charge: int = 0
    isotope: Optional[int] = None
    bonding_number: int = 0  # λ-convention
    is_aromatic: bool = False
    implicit_h: int = 0

    def __hash__(self):
        return hash(self.idx)


# =============================================================================
# P-22 MONOCYCLIC PARENT HYDRIDES / P-25 FUSED RING SYSTEMS
# =============================================================================

@dataclass
class RingInfo:
    """
    Complete information about a ring.
    Reference: P-22, P-25
    """
    atoms: Tuple[int, ...]
    size: int
    is_aromatic: bool
    heteroatoms: Dict[int, str]  # {atom_idx: element}
    element_sequence: Tuple[str, ...]
    bond_types: Tuple[int, ...] = field(default_factory=tuple)  # Bond orders

    def __hash__(self):
        return hash(self.atoms)

    @property
    def hetero_count(self) -> int:
        return len(self.heteroatoms)

    @property
    def has_nitrogen(self) -> bool:
        return 'N' in self.heteroatoms.values()

    @property
    def has_oxygen(self) -> bool:
        return 'O' in self.heteroatoms.values()

    @property
    def has_sulfur(self) -> bool:
        return 'S' in self.heteroatoms.values()

    @property
    def is_carbocyclic(self) -> bool:
        return self.hetero_count == 0

    @property
    def is_heterocyclic(self) -> bool:
        return self.hetero_count > 0


@dataclass
class FusedRingSystem:
    """
    A fused ring system (P-25).
    """
    rings: List[RingInfo]
    all_atoms: Set[int]
    fusion_bonds: List[Tuple[int, int]]  # Atoms shared between rings
    peripheral_atoms: List[int]  # Atoms on the outer edge
    interior_atoms: List[int]  # Internal fusion atoms


# =============================================================================
# P-24 SPIRO RING SYSTEMS
# =============================================================================

@dataclass
class SpiroSystem:
    """
    A spiro ring system (P-24).
    """
    rings: List[RingInfo]
    spiro_atoms: List[int]  # Atoms shared between rings (1 atom per junction)
    is_monospiro: bool
    is_polyspiro: bool
    is_branched: bool = False


# =============================================================================
# P-23 VON BAEYER SYSTEMS (BRIDGED POLYCYCLIC)
# =============================================================================

@dataclass
class VonBaeyerSystem:
    """
    A von Baeyer (bridged bicyclic/polycyclic) system (P-23).
    """
    main_ring_size: int
    bridges: List[int]  # Number of atoms in each bridge
    bridgehead_atoms: Tuple[int, int]
    total_atoms: int
    heteroatoms: Dict[int, str] = field(default_factory=dict)


# =============================================================================
# P-33 CHARACTERISTIC GROUPS (FUNCTIONAL GROUPS)
# =============================================================================

class CharacteristicGroupClass(Enum):
    """
    Classes of characteristic groups in order of seniority (P-41).
    """
    # Highest seniority
    RADICAL = 1
    ANION = 2
    CATION = 3
    ZWITTERION = 4
    # Acids
    CARBOXYLIC_ACID = 5
    SULFONIC_ACID = 6
    SULFINIC_ACID = 7
    # Derivatives
    ANHYDRIDE = 8
    ESTER = 9
    HALIDE = 10
    AMIDE = 11
    HYDRAZIDE = 12
    IMIDE = 13
    NITRILE = 14
    # Aldehydes, ketones
    ALDEHYDE = 15
    KETONE = 16
    # Alcohols, amines
    ALCOHOL = 17
    PHENOL = 18
    THIOL = 19
    AMINE = 20
    IMINE = 21
    # Others
    ETHER = 22
    PEROXIDE = 23
    # Hydrocarbons (lowest)
    HYDROCARBON = 99


@dataclass
class CharacteristicGroup:
    """
    A characteristic (functional) group.
    Reference: P-33, P-34, P-35
    """
    group_class: CharacteristicGroupClass
    name: str
    suffix: Optional[str]
    prefix: Optional[str]
    atoms: Set[int]
    attachment_atom: int
    is_principal: bool = False


# =============================================================================
# P-29, P-35 SUBSTITUENT GROUPS
# =============================================================================

@dataclass
class SubstituentGroup:
    """
    A substituent group (prefix in nomenclature).
    Reference: P-29, P-35
    """
    name: str
    atoms: Set[int]
    attachment_point: int
    locant: Union[int, str, None] = None
    is_simple: bool = True  # vs compound or complex
    multiplicity: int = 1


# =============================================================================
# P-31 UNSATURATION
# =============================================================================

@dataclass
class Unsaturation:
    """
    Double and triple bonds in the structure.
    Reference: P-31
    """
    double_bonds: List[Tuple[int, int]]
    triple_bonds: List[Tuple[int, int]]

    @property
    def total_unsaturation(self) -> int:
        return len(self.double_bonds) + 2 * len(self.triple_bonds)


# =============================================================================
# P-25.4 BRIDGES
# =============================================================================

@dataclass
class Bridge:
    """
    A bridge in a bridged fused system (P-25.4).
    """
    bridge_type: str  # 'epoxy', 'methano', 'ethano', etc.
    bridge_atoms: Set[int]
    anchor1: int
    anchor2: int
    locant1: Optional[Union[int, str]] = None
    locant2: Optional[Union[int, str]] = None


# =============================================================================
# P-90 STEREOCHEMISTRY
# =============================================================================

@dataclass
class StereocenterInfo:
    """
    Information about a stereocenter.
    Reference: P-90, P-91, P-92, P-93
    """
    atom_idx: int
    descriptor: str  # 'R', 'S', 'r', 's', 'E', 'Z', etc.
    type: str  # 'tetrahedral', 'double_bond', 'axial', etc.


# =============================================================================
# COMPLETE MOLECULAR ANALYSIS
# =============================================================================

@dataclass
class MoleculeAnalysis:
    """
    Complete analysis of a molecule for nomenclature purposes.
    """
    # Basic info
    formula: str
    total_atoms: int

    # Atoms
    atoms: List[AtomInfo]

    # Ring systems (P-22 through P-28)
    rings: List[RingInfo]
    fused_systems: List[FusedRingSystem]
    spiro_systems: List[SpiroSystem]
    von_baeyer_systems: List[VonBaeyerSystem]
    ring_assemblies: List[List[int]]  # P-28

    # Chains
    main_chain: Optional[List[int]] = None
    side_chains: List[List[int]] = field(default_factory=list)

    # Functional groups (P-33)
    characteristic_groups: List[CharacteristicGroup] = field(default_factory=list)
    principal_group: Optional[CharacteristicGroup] = None

    # Substituents (P-29, P-35)
    substituents: List[SubstituentGroup] = field(default_factory=list)

    # Unsaturation (P-31)
    unsaturation: Optional[Unsaturation] = None

    # Bridges (P-25.4)
    bridges: List[Bridge] = field(default_factory=list)

    # Stereochemistry (P-90)
    stereocenters: List[StereocenterInfo] = field(default_factory=list)

    # Parent structure type
    parent_type: str = "unknown"  # 'acyclic', 'monocyclic', 'fused', 'spiro', 'vonbaeyer', 'phane'


# =============================================================================
# NAME COMPONENTS
# =============================================================================

@dataclass
class NameComponents:
    """
    Components of an IUPAC name ready for assembly.
    """
    # Prefixes (in order)
    stereochemistry_prefix: str = ""  # (R)-, (S)-, (E)-, (Z)-, etc.
    substituent_prefixes: List[str] = field(default_factory=list)
    bridge_prefixes: List[str] = field(default_factory=list)
    fusion_prefixes: List[str] = field(default_factory=list)

    # Parent structure
    parent_name: str = ""

    # Modifications
    hydro_prefixes: List[str] = field(default_factory=list)  # P-31.2
    unsaturation_suffixes: List[str] = field(default_factory=list)  # -ene, -yne

    # Suffixes (for principal characteristic group)
    functional_suffix: str = ""

    # Locants
    locants: Dict[str, List[Union[int, str]]] = field(default_factory=dict)

    def build_name(self) -> str:
        """Assemble the complete name from components."""
        # This will be implemented in the name builder
        pass