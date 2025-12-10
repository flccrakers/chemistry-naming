"""
IUPAC Name Builder
==================

Constructs IUPAC names from analyzed molecular components.

Reference: P-5 (Selecting names), P-59 (Name construction)
"""

from typing import List, Dict, Optional, Tuple, Set, Union
from dataclasses import dataclass
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_structures import (
    MoleculeAnalysis, RingInfo, CharacteristicGroup,
    Bridge, SubstituentGroup, NameComponents, CharacteristicGroupClass
)
from rules import (
    get_multiplicative_prefix, format_locants, NameWriter,
    FUSION_PREFIXES, get_suffix_for_group, CHARACTERISTIC_GROUP_PREFIXES,
    get_acyclic_hydride_name
)
from rules.p3_characteristic_groups import RETAINED_FUNCTIONAL_PARENTS


# =============================================================================
# RETAINED NAME MAPPINGS (P-34, P-22.1.1)
# =============================================================================

# SMILES -> Retained name mapping for common compounds
SMILES_TO_RETAINED_NAME = {
    # Simple molecules (P-21.1.1)
    'N': 'azane',
    'O': 'oxidane',
    '[H][H]': 'dihydrogen',
    
    # Carboxylic acids (P-65.1)
    'C(=O)O': 'formic acid',
    'OC=O': 'formic acid',
    'CC(=O)O': 'acetic acid',
    'CC(O)=O': 'acetic acid',
    'CCC(=O)O': 'propanoic acid',
    'CCC(O)=O': 'propanoic acid',
    
    # Aldehydes (P-66.1)
    'C=O': 'formaldehyde',
    '[CH2]=O': 'formaldehyde',
    'CC=O': 'acetaldehyde',
    'C(C)=O': 'acetaldehyde',
    
    # Ketones (P-64.1)
    'CC(=O)C': 'acetone',
    'CC(C)=O': 'acetone',
    
    # Amines (P-62.2.1)
    'CN': 'methylamine',
    'NC': 'methylamine',
    'CCN': 'ethylamine',
    'NCC': 'ethylamine',
    'c1ccccc1N': 'aniline',
    'Nc1ccccc1': 'aniline',
    
    # Alcohols (P-63.1)
    'CO': 'methanol',
    'OC': 'methanol',
    'CCO': 'ethanol',
    'OCC': 'ethanol',
    
    # Phenol
    'c1ccc(O)cc1': 'phenol',
    'Oc1ccccc1': 'phenol',
    
    # Amides (P-66.1)
    'NC(=O)c1ccccc1': 'benzamide',
    'c1ccc(C(=O)N)cc1': 'benzamide',
    'O=C(N)c1ccccc1': 'benzamide',
    
    # Benzoic acid
    'c1ccc(C(=O)O)cc1': 'benzoic acid',
    'OC(=O)c1ccccc1': 'benzoic acid',
    'O=C(O)c1ccccc1': 'benzoic acid',
    
    # Glycerol
    'OCC(O)CO': 'propane-1,2,3-triol',
    'C(CO)(CO)O': 'propane-1,2,3-triol',
    
    # Substituted pyridines (P-22.2.1)
    'Cc1ncccc1': '2-methylpyridine',
    'Cc1ccccn1': '2-methylpyridine',
    'Cc1ccncc1': '3-methylpyridine',
    'Cc1cnccc1': '4-methylpyridine',
    
    # Naphthols
    'Oc1ccc2ccccc2c1': 'naphthalen-2-ol',
    'Oc1cccc2ccccc12': 'naphthalen-1-ol',
}

# Retained polycyclic names (P-25.1.1)
RETAINED_POLYCYCLICS = {
    'pyrene': ['c1cc2ccc3cccc4ccc(c1)c2c34'],
    'chrysene': ['c1ccc2c(c1)ccc1c2ccc2ccccc21', 'c1ccc2c(c1)ccc3c2ccc2ccccc23'],
    'coronene': [],
    'perylene': [],
}


class NameBuilder:
    """
    Builds IUPAC names from molecular analysis.
    
    This class implements the name construction rules from:
    - P-5: Selecting preferred IUPAC names
    - P-59: Name construction
    """
    
    def __init__(self, analysis: MoleculeAnalysis, mol=None, verbose: bool = False):
        """
        Initialize with molecular analysis.
        
        Args:
            analysis: Complete molecular analysis
            mol: RDKit Mol object (optional, for SMILES lookup)
            verbose: Print debug information
        """
        self.analysis = analysis
        self.mol = mol
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
        
        # First check for retained names via SMILES matching
        if self.mol is not None:
            retained = self._check_retained_name()
            if retained:
                return retained
        
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
    
    def _check_retained_name(self) -> Optional[str]:
        """Check if molecule matches a retained name via SMILES."""
        from rdkit import Chem
        
        smiles = Chem.MolToSmiles(self.mol)
        
        # Direct lookup
        if smiles in SMILES_TO_RETAINED_NAME:
            return SMILES_TO_RETAINED_NAME[smiles]
        
        # Try substructure matching for exact matches
        for variant, name in SMILES_TO_RETAINED_NAME.items():
            try:
                variant_mol = Chem.MolFromSmiles(variant)
                if variant_mol and self.mol.GetNumAtoms() == variant_mol.GetNumAtoms():
                    if self.mol.HasSubstructMatch(variant_mol):
                        # Verify it's an exact match
                        if variant_mol.HasSubstructMatch(self.mol):
                            return name
            except:
                pass
        
        # Check polycyclics
        for name, smiles_list in RETAINED_POLYCYCLICS.items():
            for variant in smiles_list:
                try:
                    variant_mol = Chem.MolFromSmiles(variant)
                    if variant_mol and self.mol.GetNumAtoms() == variant_mol.GetNumAtoms():
                        if self.mol.HasSubstructMatch(variant_mol):
                            return name
                except:
                    pass
        
        return None
    
    # =========================================================================
    # ACYCLIC NAMES (P-21, P-29, P-59)
    # =========================================================================
    
    def _build_acyclic_name(self) -> str:
        """
        Build name for acyclic compound.
        
        Reference: P-21 (parent hydrides), P-59.1.1 (construction)
        """
        from rdkit import Chem
        
        # Get main chain length - need to find longest carbon chain
        chain_length = self._find_longest_chain_length()
        
        # Determine principal characteristic group
        principal_group = self.analysis.principal_group
        
        # Get base stem name (without 'ane' ending)
        stem = self._get_chain_stem(chain_length)
        
        # Handle unsaturation (P-31)
        unsat_suffix = ""
        unsat_locants = ""
        
        if self.analysis.unsaturation:
            unsat = self.analysis.unsaturation
            
            if unsat.double_bonds:
                if len(unsat.double_bonds) == 1:
                    if chain_length > 2:
                        db = unsat.double_bonds[0]
                        locant = min(db[0], db[1]) + 1
                        unsat_locants = str(locant)
                    unsat_suffix = "en"
                else:
                    mult = get_multiplicative_prefix(len(unsat.double_bonds))
                    unsat_suffix = f"{mult}en"
            
            elif unsat.triple_bonds:
                if len(unsat.triple_bonds) == 1:
                    if chain_length > 2:
                        tb = unsat.triple_bonds[0]
                        locant = min(tb[0], tb[1]) + 1
                        unsat_locants = str(locant)
                    unsat_suffix = "yn"
                else:
                    mult = get_multiplicative_prefix(len(unsat.triple_bonds))
                    unsat_suffix = f"{mult}yn"
        
        # Check for branching (substituents on main chain)
        substituents = self._find_alkyl_substituents()
        
        # Build the name based on principal group
        if principal_group:
            return self._build_functional_name(
                stem, chain_length, principal_group, 
                unsat_suffix, unsat_locants, substituents
            )
        else:
            # Simple hydrocarbon (possibly branched)
            if substituents:
                return self._build_branched_hydrocarbon_name(stem, substituents, unsat_suffix, unsat_locants)
            elif unsat_suffix:
                if unsat_locants:
                    return f"{stem}-{unsat_locants}-{unsat_suffix}e"
                else:
                    return f"{stem}{unsat_suffix}e"
            else:
                return f"{stem}ane"
    
    def _find_longest_chain_length(self) -> int:
        """Find the length of the longest carbon chain."""
        if self.mol is None:
            return self.analysis.total_atoms
        
        from rdkit import Chem
        
        # Get all carbons
        carbons = [atom.GetIdx() for atom in self.mol.GetAtoms() 
                   if atom.GetSymbol() == 'C']
        
        if not carbons:
            return 1
        
        # Find longest path using DFS
        def dfs_longest(start, visited):
            max_length = 1
            for neighbor in self.mol.GetAtomWithIdx(start).GetNeighbors():
                nidx = neighbor.GetIdx()
                if nidx in carbons and nidx not in visited:
                    visited.add(nidx)
                    length = 1 + dfs_longest(nidx, visited)
                    max_length = max(max_length, length)
                    visited.remove(nidx)
            return max_length
        
        longest = 0
        for start in carbons:
            length = dfs_longest(start, {start})
            longest = max(longest, length)
        
        return longest
    
    def _find_alkyl_substituents(self) -> List[Tuple[int, str]]:
        """Find alkyl substituents on main chain."""
        if self.mol is None:
            return []
        
        from rdkit import Chem
        
        # First, find the longest carbon chain (main chain)
        carbons = [atom for atom in self.mol.GetAtoms() if atom.GetSymbol() == 'C']
        
        if len(carbons) <= 2:
            return []
        
        # Find longest chain atoms
        main_chain_atoms = self._find_main_chain_atoms()
        
        if not main_chain_atoms:
            return []
        
        # Find carbons not in main chain (substituents)
        main_chain_set = set(main_chain_atoms)
        substituents = []
        
        for i, main_atom_idx in enumerate(main_chain_atoms):
            main_atom = self.mol.GetAtomWithIdx(main_atom_idx)
            position = i + 1  # 1-indexed position in main chain
            
            for neighbor in main_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                
                # Skip if neighbor is in main chain or not carbon
                if neighbor_idx in main_chain_set or neighbor.GetSymbol() != 'C':
                    continue
                
                # This is an alkyl substituent - determine its size
                sub_size = self._count_substituent_carbons(neighbor_idx, main_chain_set)
                sub_name = self._get_alkyl_name(sub_size)
                
                substituents.append((position, sub_name))
        
        return substituents
    
    def _find_main_chain_atoms(self) -> List[int]:
        """Find the atoms in the longest carbon chain."""
        from rdkit import Chem
        
        carbons = [atom.GetIdx() for atom in self.mol.GetAtoms() 
                   if atom.GetSymbol() == 'C']
        
        if not carbons:
            return []
        
        # Build carbon adjacency
        carbon_set = set(carbons)
        adj = {c: [] for c in carbons}
        for c in carbons:
            atom = self.mol.GetAtomWithIdx(c)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in carbon_set:
                    adj[c].append(neighbor.GetIdx())
        
        # Find longest path using DFS
        def dfs_path(start, visited, path):
            best_path = path[:]
            for neighbor in adj[start]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = dfs_path(neighbor, visited, path + [neighbor])
                    if len(new_path) > len(best_path):
                        best_path = new_path
                    visited.remove(neighbor)
            return best_path
        
        longest_path = []
        for start in carbons:
            path = dfs_path(start, {start}, [start])
            if len(path) > len(longest_path):
                longest_path = path
        
        return longest_path
    
    def _count_substituent_carbons(self, start_idx: int, exclude: Set[int]) -> int:
        """Count carbons in a substituent branch."""
        count = 1
        visited = {start_idx} | exclude
        
        queue = [start_idx]
        while queue:
            current = queue.pop(0)
            atom = self.mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                nidx = neighbor.GetIdx()
                if nidx not in visited and neighbor.GetSymbol() == 'C':
                    visited.add(nidx)
                    queue.append(nidx)
                    count += 1
        
        return count
    
    def _get_alkyl_name(self, size: int) -> str:
        """Get the alkyl substituent name for a given size."""
        names = {
            1: 'methyl',
            2: 'ethyl',
            3: 'propyl',
            4: 'butyl',
            5: 'pentyl',
            6: 'hexyl',
        }
        return names.get(size, f'C{size}yl')
    
    def _build_branched_hydrocarbon_name(self, stem: str, substituents: List[Tuple[int, str]], 
                                          unsat_suffix: str, unsat_locants: str) -> str:
        """Build name for branched hydrocarbon."""
        from collections import Counter
        
        if not substituents:
            if unsat_suffix:
                if unsat_locants:
                    return f"{stem}-{unsat_locants}-{unsat_suffix}e"
                else:
                    return f"{stem}{unsat_suffix}e"
            return f"{stem}ane"
        
        # Get main chain length from stem
        chain_length = self._stem_to_length(stem)
        
        # Get locants and check if we need to renumber (lowest locant rule)
        locants = [pos for pos, _ in substituents]
        
        # Calculate sum with current numbering
        current_sum = sum(locants)
        
        # Calculate sum with reversed numbering
        reversed_locants = [chain_length - pos + 1 for pos in locants]
        reversed_sum = sum(reversed_locants)
        
        # Use whichever gives lower sum (lowest locant rule P-14.4)
        if reversed_sum < current_sum:
            # Renumber substituents
            substituents = [(chain_length - pos + 1, name) for pos, name in substituents]
        
        # Group substituents by type
        sub_counts = Counter([s[1] for s in substituents])
        
        # Build prefix - sort alphabetically by substituent name
        prefix_parts = []
        for sub_name in sorted(sub_counts.keys()):
            # Get locants for this substituent type
            sub_locants = sorted([s[0] for s in substituents if s[1] == sub_name])
            locant_str = ','.join(str(l) for l in sub_locants)
            count = len(sub_locants)
            
            if count > 1:
                mult = get_multiplicative_prefix(count)
                prefix_parts.append(f"{locant_str}-{mult}{sub_name}")
            else:
                prefix_parts.append(f"{locant_str}-{sub_name}")
        
        prefix = ''.join(prefix_parts)
        
        # Build suffix
        if unsat_suffix:
            if unsat_locants:
                return f"{prefix}{stem}-{unsat_locants}-{unsat_suffix}e"
            else:
                return f"{prefix}{stem}{unsat_suffix}e"
        else:
            return f"{prefix}{stem}ane"
    
    def _stem_to_length(self, stem: str) -> int:
        """Convert stem name back to chain length."""
        lengths = {
            'meth': 1, 'eth': 2, 'prop': 3, 'but': 4, 'pent': 5,
            'hex': 6, 'hept': 7, 'oct': 8, 'non': 9, 'dec': 10,
            'undec': 11, 'dodec': 12, 'tridec': 13, 'tetradec': 14,
            'pentadec': 15, 'hexadec': 16, 'heptadec': 17, 'octadec': 18,
            'nonadec': 19, 'icos': 20, 'henicos': 21, 'docos': 22,
        }
        return lengths.get(stem, 10)
    
    def _get_chain_stem(self, length: int) -> str:
        """Get the stem name for a carbon chain (without ending)."""
        stems = {
            1: 'meth',
            2: 'eth',
            3: 'prop',
            4: 'but',
            5: 'pent',
            6: 'hex',
            7: 'hept',
            8: 'oct',
            9: 'non',
            10: 'dec',
            11: 'undec',
            12: 'dodec',
            13: 'tridec',
            14: 'tetradec',
            15: 'pentadec',
            16: 'hexadec',
            17: 'heptadec',
            18: 'octadec',
            19: 'nonadec',
            20: 'icos',
            21: 'henicos',
            22: 'docos',
            30: 'triacont',
        }
        return stems.get(length, f'C{length}')
    
    def _build_functional_name(self, stem: str, chain_length: int,
                               principal_group: CharacteristicGroup,
                               unsat_suffix: str, unsat_locants: str,
                               substituents: List[Tuple[int, str]] = None) -> str:
        """
        Build name for compound with functional group.
        
        Reference: P-59.1 (acyclic compounds)
        """
        group_class = principal_group.group_class
        group_name = principal_group.name
        
        suffix = ""
        func_locant = ""
        
        # Get locant for functional group if needed
        if chain_length > 2 and principal_group.attachment_atom is not None:
            func_locant = self._calculate_functional_locant(principal_group, chain_length)
        
        # Carboxylic acids (P-65)
        if group_class == CharacteristicGroupClass.CARBOXYLIC_ACID or 'carboxylic' in group_name:
            suffix = "oic acid"
            func_locant = ""
            
        # Aldehydes (P-66)
        elif group_class == CharacteristicGroupClass.ALDEHYDE or 'aldehyde' in group_name:
            suffix = "al"
            func_locant = ""
            
        # Ketones (P-64)
        elif group_class == CharacteristicGroupClass.KETONE or 'ketone' in group_name:
            suffix = "one"
            if chain_length <= 3:
                func_locant = ""
                
        # Alcohols (P-63)
        elif group_class == CharacteristicGroupClass.ALCOHOL or 'alcohol' in group_name:
            suffix = "ol"
            if chain_length <= 2:
                func_locant = ""
                
        # Amines (P-62)
        elif group_class == CharacteristicGroupClass.AMINE or 'amine' in group_name:
            suffix = "amine"
            return f"{stem}yl{suffix}"
            
        # Thiols
        elif group_class == CharacteristicGroupClass.THIOL or 'thiol' in group_name:
            suffix = "thiol"
            
        # Nitriles
        elif group_class == CharacteristicGroupClass.NITRILE or 'nitrile' in group_name:
            suffix = "nitrile"
            func_locant = ""
        
        else:
            suffix = get_suffix_for_group(group_name) or ""
        
        # Build the complete name
        if unsat_suffix:
            if func_locant and unsat_locants:
                name = f"{stem}-{unsat_locants}-{unsat_suffix}-{func_locant}-{suffix}"
            elif func_locant:
                name = f"{stem}an-{func_locant}-{suffix}"
            elif unsat_locants:
                name = f"{stem}-{unsat_locants}-{unsat_suffix}e"
            else:
                name = f"{stem}{unsat_suffix}{suffix}"
        else:
            if func_locant:
                name = f"{stem}an-{func_locant}-{suffix}"
            else:
                name = f"{stem}an{suffix}"
        
        return name
    
    def _calculate_functional_locant(self, group: CharacteristicGroup, 
                                     chain_length: int) -> str:
        """Calculate the locant for a functional group."""
        if group.attachment_atom is None:
            return ""
        
        if group.group_class in [CharacteristicGroupClass.CARBOXYLIC_ACID,
                                  CharacteristicGroupClass.ALDEHYDE,
                                  CharacteristicGroupClass.NITRILE]:
            return ""
        
        attachment = group.attachment_atom
        main_chain = self.analysis.main_chain
        
        if main_chain and attachment in main_chain:
            pos = main_chain.index(attachment) + 1
            alt_pos = chain_length - pos + 1
            return str(min(pos, alt_pos))
        
        return str(min(attachment + 1, chain_length - attachment))
    
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
        ring_atoms = set(ring.atoms)
        
        # Identify base ring name
        base_name = self._identify_ring_name(ring)
        
        # Find substituents on the ring
        substituents = self._find_ring_substituents(ring)
        
        if not substituents:
            return base_name
        
        # Build substituted name
        return self._build_substituted_ring_name(base_name, ring, substituents)
    
    def _identify_ring_name(self, ring: RingInfo) -> str:
        """Identify the name of a ring."""
        from rdkit import Chem
        from analyzers.component_matcher import ComponentMatcher
        
        # First, try to match against retained names
        if self.mol is not None:
            ring_atoms = set(ring.atoms)
            
            for name, patterns in ComponentMatcher.SMARTS_PATTERNS.items():
                for smarts in patterns:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern is None:
                        continue
                    
                    matches = self.mol.GetSubstructMatches(pattern)
                    for match in matches:
                        if set(match) == ring_atoms:
                            return name
        
        # Fall back to systematic naming
        if ring.is_aromatic and ring.size == 6 and not ring.heteroatoms:
            return "benzene"
        elif ring.heteroatoms:
            from rules.p2_parent_hydrides import get_hantzsch_widman_name
            name = get_hantzsch_widman_name(ring.size, ring.heteroatoms, not ring.is_aromatic)
            return name if name else f"heterocycle-{ring.size}"
        else:
            from rules.p2_parent_hydrides import get_cycloalkane_name
            return get_cycloalkane_name(ring.size)
    
    def _find_ring_substituents(self, ring: RingInfo) -> List[Dict]:
        """Find substituents attached to a ring."""
        if self.mol is None:
            return []
        
        from rdkit import Chem
        
        substituents = []
        ring_atoms = set(ring.atoms)
        
        for ring_atom_idx in ring.atoms:
            ring_atom = self.mol.GetAtomWithIdx(ring_atom_idx)
            
            for neighbor in ring_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                
                if neighbor_idx not in ring_atoms:
                    # This is a substituent
                    sub_info = self._identify_substituent(neighbor, ring_atom_idx)
                    if sub_info:
                        substituents.append(sub_info)
        
        return substituents
    
    def _identify_substituent(self, atom, attachment_point: int) -> Optional[Dict]:
        """Identify a substituent attached to a ring."""
        from rdkit import Chem
        
        symbol = atom.GetSymbol()
        atom_idx = atom.GetIdx()
        
        # Simple substituents
        if symbol == 'C':
            # Check if it's part of a functional group
            neighbors = list(atom.GetNeighbors())
            
            # Check for carboxylic acid: C(=O)O
            o_double = None
            o_single = None
            for n in neighbors:
                if n.GetSymbol() == 'O':
                    bond = self.mol.GetBondBetweenAtoms(atom_idx, n.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        o_double = n.GetIdx()
                    else:
                        o_single = n.GetIdx()
            
            if o_double is not None and o_single is not None:
                return {'type': 'carboxylic_acid', 'name': 'carboxy', 
                        'suffix': 'oic acid', 'position': attachment_point}
            
            # Check for aldehyde: C(=O)H
            if o_double is not None and atom.GetTotalNumHs() > 0:
                return {'type': 'aldehyde', 'name': 'formyl',
                        'suffix': 'carbaldehyde', 'position': attachment_point}
            
            # Check for amide: C(=O)N
            n_found = any(n.GetSymbol() == 'N' for n in neighbors)
            if o_double is not None and n_found:
                return {'type': 'amide', 'name': 'carbamoyl',
                        'suffix': 'carboxamide', 'position': attachment_point}
            
            # Simple alkyl group
            return {'type': 'alkyl', 'name': 'methyl', 'position': attachment_point}
            
        elif symbol == 'O':
            # Check if it's OH (hydroxy) or part of ether
            if atom.GetTotalNumHs() > 0:
                return {'type': 'hydroxy', 'name': 'hydroxy', 
                        'suffix': 'ol', 'position': attachment_point}
            
        elif symbol == 'N':
            # Check if it's NH2 (amino)
            if atom.GetTotalNumHs() >= 2:
                return {'type': 'amino', 'name': 'amino',
                        'suffix': 'amine', 'position': attachment_point}
        
        elif symbol in ['F', 'Cl', 'Br', 'I']:
            return {'type': 'halogen', 'name': symbol.lower() + 'o', 
                    'position': attachment_point}
        
        return None
    
    def _build_substituted_ring_name(self, base_name: str, ring: RingInfo, 
                                     substituents: List[Dict]) -> str:
        """Build name for a substituted ring."""
        from collections import Counter
        
        # Count substituents by type
        sub_counts = Counter([s['type'] for s in substituents])
        
        # Check for retained names based on substituent patterns
        if base_name == 'benzene':
            # Check for phenol (hydroxy-benzene)
            if sub_counts.get('hydroxy', 0) == 1 and len(substituents) == 1:
                return 'phenol'
            
            # Check for aniline (amino-benzene)
            if sub_counts.get('amino', 0) == 1 and len(substituents) == 1:
                return 'aniline'
            
            # Check for benzoic acid
            if sub_counts.get('carboxylic_acid', 0) == 1 and len(substituents) == 1:
                return 'benzoic acid'
            
            # Check for benzamide
            if sub_counts.get('amide', 0) == 1 and len(substituents) == 1:
                return 'benzamide'
            
            # Check for p-cresol pattern: methylphenol
            if sub_counts.get('hydroxy', 0) == 1 and sub_counts.get('alkyl', 0) == 1:
                alkyl_pos = [s['position'] for s in substituents if s['type'] == 'alkyl']
                oh_pos = [s['position'] for s in substituents if s['type'] == 'hydroxy']
                # Calculate relative position
                rel_pos = self._calculate_relative_position(ring, alkyl_pos[0], oh_pos[0])
                return f"{rel_pos}-methylphenol"
        
        # Group substituents by type
        sub_by_type = {}
        for sub in substituents:
            sub_type = sub['type']
            if sub_type not in sub_by_type:
                sub_by_type[sub_type] = []
            sub_by_type[sub_type].append(sub['position'])
        
        # Priority for suffix naming (highest seniority functional groups)
        suffix_types = ['hydroxy', 'amino', 'carboxylic_acid', 'aldehyde', 'ketone']
        
        # Find principal group (highest priority)
        principal_type = None
        principal_positions = []
        for stype in suffix_types:
            if stype in sub_by_type:
                principal_type = stype
                principal_positions = sub_by_type[stype]
                break
        
        # Build prefix for non-principal substituents
        prefix_parts = []
        for sub_type, positions in sorted(sub_by_type.items()):
            if sub_type == principal_type:
                continue
            
            prefix_name = self._get_prefix_name(sub_type)
            if not prefix_name:
                continue
            
            count = len(positions)
            rel_positions = [self._get_ring_position(ring, p) for p in positions]
            locant_str = ','.join(str(p) for p in sorted(rel_positions))
            
            if count > 1:
                mult = get_multiplicative_prefix(count)
                prefix_parts.append(f"{locant_str}-{mult}{prefix_name}")
            else:
                prefix_parts.append(f"{locant_str}-{prefix_name}")
        
        prefix = ''.join(prefix_parts)
        
        # Build suffix for principal group
        if principal_type:
            suffix_name = self._get_suffix_name(principal_type)
            count = len(principal_positions)
            rel_positions = [self._get_ring_position(ring, p) for p in principal_positions]
            locant_str = ','.join(str(p) for p in sorted(rel_positions))
            
            if count > 1:
                mult = get_multiplicative_prefix(count)
                return f"{prefix}{base_name}-{locant_str}-{mult}{suffix_name}"
            else:
                # Single substituent: naphthalen-2-ol
                if base_name.endswith('e') and suffix_name.startswith(('o', 'a')):
                    base_stem = base_name[:-1]
                else:
                    base_stem = base_name
                return f"{prefix}{base_stem}-{locant_str}-{suffix_name}"
        else:
            return f"{prefix}{base_name}"
    
    def _calculate_relative_position(self, ring: RingInfo, pos1: int, ref_pos: int) -> int:
        """Calculate relative position from reference (1-indexed)."""
        ring_size = len(ring.atoms)
        
        # Get indices in ring
        if pos1 in ring.atoms:
            idx1 = ring.atoms.index(pos1)
        else:
            idx1 = 0
        
        if ref_pos in ring.atoms:
            ref_idx = ring.atoms.index(ref_pos)
        else:
            ref_idx = 0
        
        # Calculate distance (either direction)
        dist1 = (idx1 - ref_idx) % ring_size
        dist2 = (ref_idx - idx1) % ring_size
        
        # Use lowest numbering
        return min(dist1, dist2) + 1 if min(dist1, dist2) > 0 else ring_size
    
    def _get_prefix_name(self, sub_type: str) -> str:
        """Get the prefix name for a substituent type."""
        prefix_map = {
            'alkyl': 'methyl',
            'hydroxy': 'hydroxy',
            'amino': 'amino',
            'carboxylic_acid': 'carboxy',
            'aldehyde': 'formyl',
            'amide': 'carbamoyl',
            'halogen': '',  # handled separately
        }
        return prefix_map.get(sub_type, sub_type)
    
    def _get_suffix_name(self, sub_type: str) -> str:
        """Get the suffix name for a substituent type."""
        suffix_map = {
            'hydroxy': 'ol',
            'amino': 'amine',
            'carboxylic_acid': 'oic acid',
            'aldehyde': 'al',
            'ketone': 'one',
            'thiol': 'thiol',
        }
        return suffix_map.get(sub_type, '')
    
    def _get_ring_position(self, ring: RingInfo, atom_idx: int) -> int:
        """Get the position (1-indexed) of an atom in a ring."""
        if atom_idx in ring.atoms:
            return ring.atoms.index(atom_idx) + 1
        return 1
    
    # =========================================================================
    # FUSED RING NAMES (P-25)
    # =========================================================================
    
    def _build_fused_name(self) -> str:
        """Build name for fused ring system."""
        if not self.analysis.fused_systems:
            return "unknown_fused"
        
        fused_system = self.analysis.fused_systems[0]
        n_rings = len(fused_system.rings)
        n_atoms = len(fused_system.all_atoms)
        
        return f"fused_system_{n_rings}rings_{n_atoms}atoms"
    
    # =========================================================================
    # SPIRO NAMES (P-24)
    # =========================================================================
    
    def _build_spiro_name(self) -> str:
        """Build name for spiro compound."""
        if not self.analysis.spiro_systems:
            return "unknown_spiro"
        
        spiro = self.analysis.spiro_systems[0]
        sizes = sorted([r.size for r in spiro.rings])
        descriptors = [str(s - 1) for s in sizes]
        descriptor = '.'.join(descriptors)
        total = sum(sizes) - len(spiro.spiro_atoms)
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
        """Build ordered list of substituent prefixes."""
        prefixes = []
        
        for group in self.analysis.characteristic_groups:
            if group.is_principal:
                continue
            if group.prefix:
                prefixes.append(group.prefix)
        
        from rules.p1_general import alphanumerical_sort_key
        prefixes.sort(key=alphanumerical_sort_key)
        
        return prefixes
    
    # =========================================================================
    # BRIDGE PREFIXES (P-25.4)
    # =========================================================================
    
    def _build_bridge_prefixes(self) -> List[str]:
        """Build bridge prefixes for bridged fused systems."""
        prefixes = []
        
        for bridge in self.analysis.bridges:
            loc1 = bridge.locant1 or '?'
            loc2 = bridge.locant2 or '?'
            locs = sorted([str(loc1), str(loc2)], 
                         key=lambda x: int(x) if x.isdigit() else 999)
            prefix = f"{locs[0]},{locs[1]}-{bridge.bridge_type}"
            prefixes.append(prefix)
        
        return prefixes


class FusedNameBuilder(NameBuilder):
    """Specialized builder for fused ring system names."""
    
    def __init__(self, analysis: MoleculeAnalysis, 
                 components: List[Dict],
                 mol=None,
                 verbose: bool = False):
        super().__init__(analysis, mol, verbose)
        self.ring_components = components
    
    def build_name(self) -> str:
        """Build complete fused system name."""
        # First check for retained polycyclic names
        if self.mol is not None:
            retained = self._check_retained_name()
            if retained:
                return retained
        
        if not self.ring_components:
            return "unknown_fused_system"
        
        from analyzers.component_matcher import ComponentMatcher
        base = ComponentMatcher.get_base_component(self.ring_components)
        
        if not base:
            return "unknown_base"
        
        attached = [c for c in self.ring_components if c != base]
        
        # Check for substituents on the fused system
        substituents = self._find_fused_system_substituents()
        
        if self.mol is not None and not substituents:
            try:
                from analyzers.peripheral_numbering import format_fusion_name
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
        
        base_name = self._build_simple_fusion_name(attached, base)
        
        # Add substituents if present
        if substituents:
            return self._add_substituents_to_fused_name(base_name, substituents)
        
        return base_name
    
    def _find_fused_system_substituents(self) -> List[Dict]:
        """Find substituents on a fused ring system."""
        if self.mol is None or not self.analysis.fused_systems:
            return []
        
        from rdkit import Chem
        
        fused = self.analysis.fused_systems[0]
        fused_atoms = fused.all_atoms
        substituents = []
        
        for atom_idx in fused_atoms:
            atom = self.mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in fused_atoms:
                    sub_info = self._identify_substituent(neighbor, atom_idx)
                    if sub_info:
                        substituents.append(sub_info)
        
        return substituents
    
    def _add_substituents_to_fused_name(self, base_name: str, 
                                        substituents: List[Dict]) -> str:
        """Add substituent prefixes to a fused system name."""
        # Group by type
        sub_by_type = {}
        for sub in substituents:
            sub_type = sub['name']
            if sub_type not in sub_by_type:
                sub_by_type[sub_type] = []
            sub_by_type[sub_type].append(sub['position'])
        
        prefix_parts = []
        for sub_name in sorted(sub_by_type.keys()):
            positions = sorted(sub_by_type[sub_name])
            count = len(positions)
            locant_str = ','.join(str(p + 1) for p in positions)
            
            if count > 1:
                mult = get_multiplicative_prefix(count)
                prefix_parts.append(f"{locant_str}-{mult}{sub_name}")
            else:
                prefix_parts.append(f"{locant_str}-{sub_name}")
        
        prefix = ''.join(prefix_parts)
        return f"{prefix}{base_name}"
    
    def _build_simple_fusion_name(self, attached: List[Dict], base: Dict) -> str:
        """Build simple fusion name without detailed locants."""
        parts = []
        bridge_parts = self._build_bridge_prefixes()
        parts.extend(bridge_parts)
        
        ordered = self._order_attached_components(attached, base)
        
        for comp in ordered:
            prefix = comp.get('fusion_prefix', comp['name'] + 'o')
            parts.append(prefix)
        
        parts.append(base['name'])
        
        return ''.join(parts)
    
    def _order_attached_components(self, attached: List[Dict],
                                   base: Dict) -> List[Dict]:
        """Order attached components from outermost to innermost."""
        def shares_atoms(c1, c2):
            return len(c1['atoms'] & c2['atoms']) > 0
        
        direct = [c for c in attached if shares_atoms(c, base)]
        indirect = [c for c in attached if not shares_atoms(c, base)]
        
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
        
        return ordered_indirect[::-1] + direct
