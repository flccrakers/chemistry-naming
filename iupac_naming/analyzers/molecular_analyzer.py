"""
Molecular Analyzer
==================

Analyzes molecular structures using RDKit to extract information
needed for IUPAC nomenclature.

This module provides the bridge between RDKit molecules and the
IUPAC nomenclature rules.
"""

from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict
import sys
import os

# RDKit imports
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, Descriptors

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_structures import (
    AtomInfo, RingInfo, FusedRingSystem, SpiroSystem, VonBaeyerSystem,
    CharacteristicGroup, CharacteristicGroupClass, SubstituentGroup,
    Unsaturation, Bridge, StereocenterInfo, MoleculeAnalysis
)


class MolecularAnalyzer:
    """
    Analyzes a molecule and extracts structural information
    for IUPAC nomenclature.
    """
    
    def __init__(self, mol: Chem.Mol):
        """
        Initialize with an RDKit molecule.
        
        Args:
            mol: RDKit Mol object
        """
        self.mol = mol
        self._rings: Optional[List[RingInfo]] = None
        self._atoms: Optional[List[AtomInfo]] = None
        
    # =========================================================================
    # ATOM ANALYSIS
    # =========================================================================
    
    def get_atoms(self) -> List[AtomInfo]:
        """Extract information about all atoms."""
        if self._atoms is not None:
            return self._atoms
        
        self._atoms = []
        for atom in self.mol.GetAtoms():
            info = AtomInfo(
                idx=atom.GetIdx(),
                symbol=atom.GetSymbol(),
                charge=atom.GetFormalCharge(),
                isotope=atom.GetIsotope() if atom.GetIsotope() != 0 else None,
                is_aromatic=atom.GetIsAromatic(),
                implicit_h=atom.GetTotalNumHs()
            )
            self._atoms.append(info)
        
        return self._atoms
    
    def get_atom_symbol(self, idx: int) -> str:
        """Get symbol for atom at index."""
        return self.mol.GetAtomWithIdx(idx).GetSymbol()
    
    def get_neighbors(self, idx: int) -> List[int]:
        """Get indices of neighboring atoms."""
        atom = self.mol.GetAtomWithIdx(idx)
        return [n.GetIdx() for n in atom.GetNeighbors()]
    
    # =========================================================================
    # RING ANALYSIS (P-22, P-25)
    # =========================================================================
    
    def get_rings(self) -> List[RingInfo]:
        """
        Extract all rings from the molecule.
        
        Reference: P-22 (monocyclic), P-25 (fused)
        """
        if self._rings is not None:
            return self._rings
        
        ring_info = self.mol.GetRingInfo()
        self._rings = []
        
        for ring_atoms in ring_info.AtomRings():
            atoms_tuple = tuple(ring_atoms)
            
            # Get element sequence
            elements = tuple(
                self.mol.GetAtomWithIdx(idx).GetSymbol() 
                for idx in ring_atoms
            )
            
            # Check aromaticity
            is_aromatic = all(
                self.mol.GetAtomWithIdx(idx).GetIsAromatic() 
                for idx in ring_atoms
            )
            
            # Find heteroatoms
            heteroatoms = {}
            for idx in ring_atoms:
                symbol = self.mol.GetAtomWithIdx(idx).GetSymbol()
                if symbol not in ['C', 'H']:
                    heteroatoms[idx] = symbol
            
            # Get bond types in ring
            bond_types = []
            for i in range(len(ring_atoms)):
                a1 = ring_atoms[i]
                a2 = ring_atoms[(i + 1) % len(ring_atoms)]
                bond = self.mol.GetBondBetweenAtoms(a1, a2)
                if bond:
                    bond_types.append(int(bond.GetBondType()))
            
            ring = RingInfo(
                atoms=atoms_tuple,
                size=len(ring_atoms),
                is_aromatic=is_aromatic,
                heteroatoms=heteroatoms,
                element_sequence=elements,
                bond_types=tuple(bond_types)
            )
            self._rings.append(ring)
        
        return self._rings
    
    def get_ring_system_atoms(self) -> Set[int]:
        """Get all atoms that are part of any ring."""
        rings = self.get_rings()
        all_atoms = set()
        for ring in rings:
            all_atoms.update(ring.atoms)
        return all_atoms
    
    def find_fused_rings(self) -> Dict[int, List[int]]:
        """
        Find which rings share atoms (are fused).
        
        Reference: P-25.3
        """
        rings = self.get_rings()
        fusion_map = defaultdict(list)
        
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                shared = set(rings[i].atoms) & set(rings[j].atoms)
                if len(shared) >= 2:  # Ortho-fusion
                    fusion_map[i].append(j)
                    fusion_map[j].append(i)
        
        return fusion_map
    
    def get_fused_ring_systems(self) -> List[FusedRingSystem]:
        """
        Identify fused ring systems.
        
        Reference: P-25
        """
        rings = self.get_rings()
        fusion_map = self.find_fused_rings()
        
        if not fusion_map:
            return []
        
        # Find connected components of fused rings
        visited = set()
        systems = []
        
        for start_ring in range(len(rings)):
            if start_ring in visited:
                continue
            if start_ring not in fusion_map:
                continue
            
            # BFS to find all connected fused rings
            component_rings = []
            queue = [start_ring]
            
            while queue:
                ring_idx = queue.pop(0)
                if ring_idx in visited:
                    continue
                visited.add(ring_idx)
                component_rings.append(rings[ring_idx])
                
                for neighbor in fusion_map.get(ring_idx, []):
                    if neighbor not in visited:
                        queue.append(neighbor)
            
            if len(component_rings) > 1:
                # Build fused system
                all_atoms = set()
                for ring in component_rings:
                    all_atoms.update(ring.atoms)
                
                # Find fusion bonds (shared atoms)
                fusion_bonds = []
                for i, r1 in enumerate(component_rings):
                    for j, r2 in enumerate(component_rings[i+1:], i+1):
                        shared = list(set(r1.atoms) & set(r2.atoms))
                        if len(shared) >= 2:
                            fusion_bonds.append(tuple(sorted(shared[:2])))
                
                system = FusedRingSystem(
                    rings=component_rings,
                    all_atoms=all_atoms,
                    fusion_bonds=fusion_bonds,
                    peripheral_atoms=[],  # TODO: calculate
                    interior_atoms=[]  # TODO: calculate
                )
                systems.append(system)
        
        return systems
    
    def find_spiro_systems(self) -> List[SpiroSystem]:
        """
        Identify spiro ring systems.
        
        Reference: P-24
        """
        rings = self.get_rings()
        spiro_atoms = []
        
        # Find atoms shared by exactly 2 rings (single atom = spiro)
        atom_ring_count = defaultdict(list)
        for i, ring in enumerate(rings):
            for atom in ring.atoms:
                atom_ring_count[atom].append(i)
        
        for atom, ring_indices in atom_ring_count.items():
            if len(ring_indices) == 2:
                # Check it's the only shared atom between these rings
                r1, r2 = ring_indices
                shared = set(rings[r1].atoms) & set(rings[r2].atoms)
                if len(shared) == 1:
                    spiro_atoms.append(atom)
        
        if not spiro_atoms:
            return []
        
        # Group into spiro systems
        # For now, return single system with all spiro atoms
        spiro_system = SpiroSystem(
            rings=rings,
            spiro_atoms=spiro_atoms,
            is_monospiro=len(spiro_atoms) == 1,
            is_polyspiro=len(spiro_atoms) > 1,
            is_branched=False  # TODO: detect branching
        )
        
        return [spiro_system]
    
    # =========================================================================
    # FUNCTIONAL GROUP ANALYSIS (P-33, P-34, P-35)
    # =========================================================================
    
    def find_characteristic_groups(self) -> List[CharacteristicGroup]:
        """
        Find characteristic (functional) groups in the molecule.
        
        Reference: P-33, P-34, P-35
        """
        groups = []
        
        # SMARTS patterns for functional groups
        patterns = {
            # Acids
            'carboxylic_acid': '[CX3](=O)[OX2H1]',
            'sulfonic_acid': '[SX4](=O)(=O)[OX2H1]',
            
            # Carbonyl compounds
            'aldehyde': '[CX3H1](=O)[#6,H]',
            'ketone': '[#6][CX3](=O)[#6]',
            
            # Alcohols and thiols
            'alcohol': '[OX2H][CX4]',
            'phenol': '[OX2H][cX3]',
            'thiol': '[SX2H]',
            
            # Amines
            'primary_amine': '[NX3H2][CX4]',
            'secondary_amine': '[NX3H1]([CX4])[CX4]',
            'tertiary_amine': '[NX3]([CX4])([CX4])[CX4]',
            
            # Nitriles
            'nitrile': '[NX1]#[CX2]',
            
            # Nitro
            'nitro': '[NX3](=O)=O',
            
            # Halides
            'fluoro': '[FX1]',
            'chloro': '[ClX1]',
            'bromo': '[BrX1]',
            'iodo': '[IX1]',
            
            # Ethers
            'ether': '[OX2]([CX4])[CX4]',
            
            # Esters
            'ester': '[CX3](=O)[OX2][CX4]',
            
            # Amides
            'amide': '[CX3](=O)[NX3]',
        }
        
        class_mapping = {
            'carboxylic_acid': CharacteristicGroupClass.CARBOXYLIC_ACID,
            'sulfonic_acid': CharacteristicGroupClass.SULFONIC_ACID,
            'aldehyde': CharacteristicGroupClass.ALDEHYDE,
            'ketone': CharacteristicGroupClass.KETONE,
            'alcohol': CharacteristicGroupClass.ALCOHOL,
            'phenol': CharacteristicGroupClass.PHENOL,
            'thiol': CharacteristicGroupClass.THIOL,
            'primary_amine': CharacteristicGroupClass.AMINE,
            'secondary_amine': CharacteristicGroupClass.AMINE,
            'tertiary_amine': CharacteristicGroupClass.AMINE,
            'nitrile': CharacteristicGroupClass.NITRILE,
            'ether': CharacteristicGroupClass.ETHER,
            'ester': CharacteristicGroupClass.ESTER,
            'amide': CharacteristicGroupClass.AMIDE,
        }
        
        suffix_mapping = {
            'carboxylic_acid': 'oic acid',
            'aldehyde': 'al',
            'ketone': 'one',
            'alcohol': 'ol',
            'thiol': 'thiol',
            'primary_amine': 'amine',
            'nitrile': 'nitrile',
            'amide': 'amide',
        }
        
        prefix_mapping = {
            'carboxylic_acid': 'carboxy',
            'aldehyde': 'formyl',
            'ketone': 'oxo',
            'alcohol': 'hydroxy',
            'phenol': 'hydroxy',
            'thiol': 'sulfanyl',
            'primary_amine': 'amino',
            'secondary_amine': 'amino',
            'tertiary_amine': 'amino',
            'nitrile': 'cyano',
            'nitro': 'nitro',
            'fluoro': 'fluoro',
            'chloro': 'chloro',
            'bromo': 'bromo',
            'iodo': 'iodo',
            'ether': None,
            'ester': None,
            'amide': 'carbamoyl',
        }
        
        for name, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            
            matches = self.mol.GetSubstructMatches(pattern)
            for match in matches:
                # Determine the correct attachment atom (should be the carbon)
                # For most patterns, the carbon is a specific position
                attachment = match[0]
                
                if name == 'alcohol':
                    # Pattern [OX2H][CX4] - carbon is at position 1
                    attachment = match[1] if len(match) > 1 else match[0]
                elif name == 'phenol':
                    # Pattern [OX2H][cX3] - carbon is at position 1
                    attachment = match[1] if len(match) > 1 else match[0]
                elif name == 'thiol':
                    # Pattern [SX2H] - find attached carbon
                    s_atom = self.mol.GetAtomWithIdx(match[0])
                    for neighbor in s_atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C':
                            attachment = neighbor.GetIdx()
                            break
                elif name in ['primary_amine', 'secondary_amine', 'tertiary_amine']:
                    # Pattern starts with N, carbon is at position 1
                    attachment = match[1] if len(match) > 1 else match[0]
                elif name == 'carboxylic_acid':
                    # Pattern [CX3](=O)[OX2H1] - C is at position 0
                    attachment = match[0]
                elif name == 'aldehyde':
                    # Pattern [CX3H1](=O) - C is at position 0
                    attachment = match[0]
                elif name == 'ketone':
                    # Pattern [#6][CX3](=O)[#6] - carbonyl C is at position 1
                    attachment = match[1] if len(match) > 1 else match[0]
                elif name == 'nitrile':
                    # Pattern [NX1]#[CX2] - C is at position 1
                    attachment = match[1] if len(match) > 1 else match[0]
                
                group = CharacteristicGroup(
                    group_class=class_mapping.get(name, CharacteristicGroupClass.HYDROCARBON),
                    name=name,
                    suffix=suffix_mapping.get(name),
                    prefix=prefix_mapping.get(name),
                    atoms=set(match),
                    attachment_atom=attachment,
                    is_principal=False
                )
                groups.append(group)
        
        return groups
    
    # =========================================================================
    # UNSATURATION ANALYSIS (P-31)
    # =========================================================================
    
    def find_unsaturation(self) -> Unsaturation:
        """
        Find double and triple bonds (excluding C=O in functional groups).
        
        Reference: P-31
        """
        double_bonds = []
        triple_bonds = []
        
        # First identify carbonyl carbons (C=O in acids, aldehydes, ketones, etc.)
        # These should NOT be counted as unsaturation
        carbonyl_bonds = set()
        
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        bond = self.mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            # This is a C=O carbonyl bond - exclude it
                            carbonyl_bonds.add((atom.GetIdx(), neighbor.GetIdx()))
                            carbonyl_bonds.add((neighbor.GetIdx(), atom.GetIdx()))
        
        for bond in self.mol.GetBonds():
            bond_type = bond.GetBondType()
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            
            # Skip aromatic bonds (handled separately)
            if bond.GetIsAromatic():
                continue
            
            # Skip carbonyl bonds (C=O in functional groups)
            if (a1, a2) in carbonyl_bonds:
                continue
            
            if bond_type == Chem.BondType.DOUBLE:
                double_bonds.append((a1, a2))
            elif bond_type == Chem.BondType.TRIPLE:
                triple_bonds.append((a1, a2))
        
        return Unsaturation(
            double_bonds=double_bonds,
            triple_bonds=triple_bonds
        )
    
    # =========================================================================
    # BRIDGE ANALYSIS (P-25.4)
    # =========================================================================
    
    def find_bridges(self) -> List[Bridge]:
        """
        Find bridges in fused ring systems.
        
        Reference: P-25.4
        """
        bridges = []
        rings = self.get_rings()
        ring_atoms = self.get_ring_system_atoms()
        
        # Find epoxy bridges
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() != 'O':
                continue
            
            o_idx = atom.GetIdx()
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            
            if len(neighbors) != 2:
                continue
            
            c1, c2 = neighbors
            
            # Both must be carbons
            if self.get_atom_symbol(c1) != 'C':
                continue
            if self.get_atom_symbol(c2) != 'C':
                continue
            
            # Not directly bonded
            direct_bond = self.mol.GetBondBetweenAtoms(c1, c2)
            if direct_bond is not None:
                continue
            
            # Find path without oxygen
            from collections import deque
            queue = deque([(c1, 0)])
            visited = {c1, o_idx}
            found_path = -1
            
            while queue:
                node, dist = queue.popleft()
                if node == c2:
                    found_path = dist
                    break
                
                for neighbor in self.get_neighbors(node):
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append((neighbor, dist + 1))
            
            if found_path >= 3:
                bridge = Bridge(
                    bridge_type='epoxy',
                    bridge_atoms={o_idx},
                    anchor1=c1,
                    anchor2=c2
                )
                bridges.append(bridge)
        
        return bridges
    
    # =========================================================================
    # STEREOCHEMISTRY ANALYSIS (P-90)
    # =========================================================================
    
    def find_stereocenters(self) -> List[StereocenterInfo]:
        """
        Find stereocenters in the molecule.
        
        Reference: P-90, P-91, P-92, P-93
        """
        stereocenters = []
        
        # Find chiral centers
        chiral_centers = Chem.FindMolChiralCenters(self.mol, includeUnassigned=True)
        
        for atom_idx, descriptor in chiral_centers:
            info = StereocenterInfo(
                atom_idx=atom_idx,
                descriptor=descriptor,
                type='tetrahedral'
            )
            stereocenters.append(info)
        
        # Find E/Z double bonds
        for bond in self.mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo in [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]:
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    descriptor = 'E' if stereo == Chem.BondStereo.STEREOE else 'Z'
                    info = StereocenterInfo(
                        atom_idx=a1,  # Use first atom as reference
                        descriptor=descriptor,
                        type='double_bond'
                    )
                    stereocenters.append(info)
        
        return stereocenters
    
    # =========================================================================
    # CHAIN ANALYSIS
    # =========================================================================
    
    def find_longest_chain(self) -> List[int]:
        """
        Find the longest carbon chain.
        
        Reference: P-44.3
        """
        # Get all carbon atoms not in rings
        ring_atoms = self.get_ring_system_atoms()
        chain_carbons = []
        
        for atom in self.mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetIdx() not in ring_atoms:
                chain_carbons.append(atom.GetIdx())
        
        if not chain_carbons:
            return []
        
        # Find longest path using DFS
        def dfs(start, visited):
            max_path = [start]
            for neighbor in self.get_neighbors(start):
                if neighbor in chain_carbons and neighbor not in visited:
                    path = dfs(neighbor, visited | {start})
                    if len(path) + 1 > len(max_path):
                        max_path = [start] + path
            return max_path
        
        longest = []
        for start in chain_carbons:
            path = dfs(start, set())
            if len(path) > len(longest):
                longest = path
        
        return longest
    
    # =========================================================================
    # COMPLETE ANALYSIS
    # =========================================================================
    
    def analyze(self) -> MoleculeAnalysis:
        """
        Perform complete analysis of the molecule.
        
        Returns:
            MoleculeAnalysis with all structural information
        """
        # Basic info
        formula = rdMolDescriptors.CalcMolFormula(self.mol)
        
        # Get all components
        atoms = self.get_atoms()
        rings = self.get_rings()
        fused_systems = self.get_fused_ring_systems()
        spiro_systems = self.find_spiro_systems()
        characteristic_groups = self.find_characteristic_groups()
        unsaturation = self.find_unsaturation()
        bridges = self.find_bridges()
        stereocenters = self.find_stereocenters()
        main_chain = self.find_longest_chain()
        
        # Determine principal characteristic group
        principal_group = None
        if characteristic_groups:
            # Sort by seniority and select highest
            from rules.p4_name_construction import get_class_seniority
            characteristic_groups.sort(
                key=lambda g: get_class_seniority(g.name),
                reverse=True
            )
            principal_group = characteristic_groups[0]
            principal_group.is_principal = True
        
        # Determine parent type
        # Special case: a system can be both "fused" (within components) and "spiro" (between components)
        # For example, trispiro with cyclopenta[b]pyran has internal fusion but is globally spiro
        
        has_spiro_junctions = spiro_systems and len(spiro_systems) > 0 and any(
            s.is_polyspiro or len(s.spiro_atoms) >= 1 for s in spiro_systems
        )
        
        # Count true spiro connections (single atoms connecting distinct ring groups)
        ring_info = self.mol.GetRingInfo()
        all_rings = list(ring_info.AtomRings())
        
        if has_spiro_junctions and len(all_rings) >= 3:
            # Check if this is primarily a spiro system
            # Count atoms in exactly 2 rings that are NOT adjacent to other shared atoms
            atom_ring_count = {}
            for ring in all_rings:
                for atom_idx in ring:
                    atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1
            
            # True spiro atoms: in 2+ rings, and no adjacent atom is also in 2+ rings
            true_spiro_count = 0
            fused_shared_count = 0
            
            for atom_idx, count in atom_ring_count.items():
                if count >= 2:
                    atom = self.mol.GetAtomWithIdx(atom_idx)
                    has_shared_neighbor = False
                    for neighbor in atom.GetNeighbors():
                        n_idx = neighbor.GetIdx()
                        if atom_ring_count.get(n_idx, 0) >= 2:
                            has_shared_neighbor = True
                            break
                    
                    if has_shared_neighbor:
                        fused_shared_count += 1
                    else:
                        true_spiro_count += 1
            
            # If we have true spiro junctions connecting groups, it's primarily spiro
            if true_spiro_count >= 2:
                parent_type = 'spiro'
            elif fused_systems:
                parent_type = 'fused'
            else:
                parent_type = 'spiro'
        elif fused_systems:
            parent_type = 'fused'
        elif spiro_systems:
            parent_type = 'spiro'
        elif rings:
            parent_type = 'monocyclic' if len(rings) == 1 else 'ring_assembly'
        elif main_chain:
            parent_type = 'acyclic'
        else:
            parent_type = 'unknown'
        
        return MoleculeAnalysis(
            formula=formula,
            total_atoms=self.mol.GetNumAtoms(),
            atoms=atoms,
            rings=rings,
            fused_systems=fused_systems,
            spiro_systems=spiro_systems,
            von_baeyer_systems=[],  # TODO
            ring_assemblies=[],  # TODO
            main_chain=main_chain,
            characteristic_groups=characteristic_groups,
            principal_group=principal_group,
            unsaturation=unsaturation,
            bridges=bridges,
            stereocenters=stereocenters,
            parent_type=parent_type
        )
