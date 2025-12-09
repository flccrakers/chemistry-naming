"""
IUPAC Name Builder
==================

Constructs IUPAC names from analyzed molecular components.

Reference: P-5 (Selecting names), P-59 (Name construction)
"""

from typing import List, Dict
import sys
import os


from ..analyzers import format_fusion_name, ComponentMatcher
from ..data_structures import NameComponents, MoleculeAnalysis
from ..rules import get_acyclic_hydride_name, get_suffix_for_group, get_hantzsch_widman_name, \
    get_cycloalkane_name, alphanumerical_sort_key

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class NameBuilder:
    """
    Builds IUPAC names from molecular analysis.

    This class implements the name construction rules from:
    - P-5: Selecting preferred IUPAC names
    - P-59: Name construction
    """

    def __init__(self, analysis: MoleculeAnalysis, verbose: bool = False):
        """
        Initialize with molecular analysis.

        Args:
            analysis: Complete molecular analysis
            verbose: Print debug information
        """
        self.analysis = analysis
        self.verbose = verbose
        self.components = NameComponents()

    def _log(self, msg: str):
        if self.verbose:
            print(f"  [NameBuilder] {msg}")

    # =========================================================================
    # MAIN NAME BUILDING
    # =========================================================================

    def build_name(self) -> str:
        """
        Build the complete IUPAC name.

        Returns:
            IUPAC name string
        """
        self._log(f"Building name for {self.analysis.parent_type} structure")

        # Route to appropriate builder based on structure type
        if self.analysis.parent_type == 'fused':
            return self._build_fused_name()
        elif self.analysis.parent_type == 'spiro':
            return self._build_spiro_name()
        elif self.analysis.parent_type == 'monocyclic':
            return self._build_monocyclic_name()
        elif self.analysis.parent_type == 'acyclic':
            return self._build_acyclic_name()
        else:
            return self._build_generic_name()

    # =========================================================================
    # ACYCLIC NAMES (P-21, P-29)
    # =========================================================================

    def _build_acyclic_name(self) -> str:
        """
        Build name for acyclic compound.

        Reference: P-21 (parent hydrides), P-59.1.1 (construction)
        """
        parts = []

        # Get main chain length
        chain = self.analysis.main_chain
        if not chain:
            return "unknown"

        chain_length = len(chain)

        # Get base name
        base_name = get_acyclic_hydride_name('C', chain_length)

        # Handle unsaturation (P-31)
        if self.analysis.unsaturation:
            unsat = self.analysis.unsaturation

            # Remove 'ane' ending
            if base_name.endswith('ane'):
                stem = base_name[:-3]
            else:
                stem = base_name

            # Add unsaturation
            if unsat.triple_bonds:
                stem += 'yn'
            elif unsat.double_bonds:
                stem += 'en'
            else:
                stem += 'an'

            base_name = stem + 'e'

        # Handle principal characteristic group (suffix)
        if self.analysis.principal_group:
            pg = self.analysis.principal_group
            suffix = get_suffix_for_group(pg.name)

            if suffix:
                # Remove final 'e' before suffix starting with vowel
                if base_name.endswith('e') and suffix[0] in 'aeiou':
                    base_name = base_name[:-1]
                base_name += suffix

        # Build substituent prefixes
        substituent_parts = self._build_substituent_prefixes()

        # Assemble name
        name = ''.join(substituent_parts) + base_name

        return name

    # =========================================================================
    # MONOCYCLIC NAMES (P-22)
    # =========================================================================

    def _build_monocyclic_name(self) -> str:
        """
        Build name for monocyclic compound.

        Reference: P-22
        """
        if not self.analysis.rings:
            return "unknown"

        ring = self.analysis.rings[0]

        # Try to identify the ring
        # For now, use simple naming
        if ring.is_aromatic and ring.size == 6 and not ring.heteroatoms:
            base_name = "benzene"
        elif ring.is_heterocyclic:
            # Use Hantzsch-Widman or retained name

            base_name = get_hantzsch_widman_name(
                ring.size, ring.heteroatoms, not ring.is_aromatic
            )
            if not base_name:
                base_name = f"heterocycle-{ring.size}"
        else:
            # Cycloalkane

            base_name = get_cycloalkane_name(ring.size)

        # Add substituents
        substituent_parts = self._build_substituent_prefixes()

        return ''.join(substituent_parts) + base_name

    # =========================================================================
    # FUSED RING NAMES (P-25)
    # =========================================================================

    def _build_fused_name(self) -> str:
        """
        Build name for fused ring system.

        Reference: P-25.3 (constructing fusion names)
        """
        if not self.analysis.fused_systems:
            return "unknown_fused"

        fused_system = self.analysis.fused_systems[0]

        # This is a placeholder - full implementation requires:
        # 1. Component identification
        # 2. Base component selection
        # 3. Fusion locant calculation
        # 4. Peripheral numbering

        # For now, return placeholder
        n_rings = len(fused_system.rings)
        n_atoms = len(fused_system.all_atoms)

        return f"fused_system_{n_rings}rings_{n_atoms}atoms"

    # =========================================================================
    # SPIRO NAMES (P-24)
    # =========================================================================

    def _build_spiro_name(self) -> str:
        """
        Build name for spiro compound.

        Reference: P-24
        """
        if not self.analysis.spiro_systems:
            return "unknown_spiro"

        spiro = self.analysis.spiro_systems[0]

        # Get ring sizes
        sizes = sorted([r.size for r in spiro.rings])

        # Build descriptor
        # spiro[a.b]parent where a,b are ring sizes minus 1
        descriptors = [str(s - 1) for s in sizes]
        descriptor = '.'.join(descriptors)

        # Get total atoms
        total = sum(sizes) - len(spiro.spiro_atoms)

        # Get base name
        base = get_acyclic_hydride_name('C', total)

        return f"spiro[{descriptor}]{base}"

    # =========================================================================
    # GENERIC NAMES
    # =========================================================================

    def _build_generic_name(self) -> str:
        """Build a generic name when specific type is unknown."""
        return f"compound_C{self.analysis.total_atoms}"

    # =========================================================================
    # SUBSTITUENT PREFIXES (P-29, P-35)
    # =========================================================================

    def _build_substituent_prefixes(self) -> List[str]:
        """
        Build ordered list of substituent prefixes.

        Reference: P-29, P-35, P-14.5 (alphanumerical order)
        """
        prefixes = []

        # Get non-principal characteristic groups
        for group in self.analysis.characteristic_groups:
            if group.is_principal:
                continue

            if group.prefix:
                prefixes.append(group.prefix)

        # Sort alphabetically (P-14.5)

        prefixes.sort(key=alphanumerical_sort_key)

        return prefixes

    # =========================================================================
    # BRIDGE PREFIXES (P-25.4)
    # =========================================================================

    def _build_bridge_prefixes(self) -> List[str]:
        """
        Build bridge prefixes for bridged fused systems.

        Reference: P-25.4
        """
        prefixes = []

        for bridge in self.analysis.bridges:
            loc1 = bridge.locant1 or '?'
            loc2 = bridge.locant2 or '?'

            # Sort locants
            locs = sorted([str(loc1), str(loc2)],
                          key=lambda x: int(x) if x.isdigit() else 999)

            prefix = f"{locs[0]},{locs[1]}-{bridge.bridge_type}"
            prefixes.append(prefix)

        return prefixes

    # =========================================================================
    # FUSION PREFIXES (P-25.3)
    # =========================================================================

    def _build_fusion_prefixes(self, components: List[Dict],
                               base: Dict) -> List[str]:
        """
        Build fusion prefixes with locants.

        Reference: P-25.3.2, P-25.3.4

        Format: prefix[locants]
        Example: pyrido[1'',2'':1',2']
        """
        prefixes = []

        for comp in components:
            if comp == base:
                continue

            fusion_prefix = comp.get('fusion_prefix', comp['name'] + 'o')

            # Calculate fusion locants
            # This requires peripheral numbering - placeholder for now
            locants = self._calculate_fusion_locants(comp, base)

            if locants:
                prefixes.append(f"{fusion_prefix}{locants}")
            else:
                prefixes.append(fusion_prefix)

        return prefixes

    def _calculate_fusion_locants(self, attached: Dict, base: Dict) -> str:
        """
        Calculate fusion locant string.

        Reference: P-25.3.2

        Returns string like [4',5':5,6] or [2,3-b]
        """
        # Find shared atoms
        shared = attached['atoms'] & base['atoms']

        if not shared:
            return ""

        # This is a placeholder - real implementation needs
        # peripheral numbering of both components
        return ""


class FusedNameBuilder(NameBuilder):
    """
    Specialized builder for fused ring system names.

    Reference: P-25
    """

    def __init__(self, analysis: MoleculeAnalysis,
                 components: List[Dict],
                 mol=None,
                 verbose: bool = False):
        super().__init__(analysis, verbose)
        self.ring_components = components
        self.mol = mol

    def build_name(self) -> str:
        """Build complete fused system name."""
        if not self.ring_components:
            return "unknown_fused_system"

        # Select base component (P-25.3.4.2.1)

        base = ComponentMatcher.get_base_component(self.ring_components)

        if not base:
            return "unknown_base"

        # Get attached components
        attached = [c for c in self.ring_components if c != base]

        # If we have the mol object, use full fusion naming
        if self.mol is not None:
            try:

                # Prepare bridges
                bridges = []
                for bridge in self.analysis.bridges:
                    bridges.append({
                        'bridge_type': bridge.bridge_type,
                        'locant1': bridge.locant1 or '?',
                        'locant2': bridge.locant2 or '?',
                    })

                return format_fusion_name(attached, base, self.mol, bridges)
            except Exception as e:
                self._log(f"Full fusion naming failed: {e}")
                # Fall through to simple naming

        # Simple naming without full locants
        return self._build_simple_fusion_name(attached, base)

    def _build_simple_fusion_name(self, attached: List[Dict], base: Dict) -> str:
        """Build simple fusion name without detailed locants."""
        parts = []

        # Add bridges
        bridge_parts = self._build_bridge_prefixes()
        parts.extend(bridge_parts)

        # Order attached components (outermost first)
        ordered = self._order_attached_components(attached, base)

        # Add fusion prefixes
        for comp in ordered:
            prefix = comp.get('fusion_prefix', comp['name'] + 'o')
            parts.append(prefix)

        # Add base
        parts.append(base['name'])

        return ''.join(parts)

    def _order_attached_components(self, attached: List[Dict],
                                   base: Dict) -> List[Dict]:
        """
        Order attached components from outermost to innermost.

        Reference: P-25.3.4.2.3
        """

        # Calculate distance from base for each component
        def shares_atoms(c1, c2):
            return len(c1['atoms'] & c2['atoms']) > 0

        # Separate into direct and indirect attachments
        direct = [c for c in attached if shares_atoms(c, base)]
        indirect = [c for c in attached if not shares_atoms(c, base)]

        # Order indirect by attachment chain
        ordered_indirect = []
        remaining = indirect.copy()
        current_level = direct

        while remaining:
            next_level = []
            for comp in remaining:
                for attached_comp in current_level:
                    if shares_atoms(comp, attached_comp):
                        next_level.append(comp)
                        break

            if next_level:
                ordered_indirect.extend(next_level)
                for c in next_level:
                    if c in remaining:
                        remaining.remove(c)
                current_level = next_level
            else:
                break

        # Final order: outermost first
        return ordered_indirect[::-1] + direct
