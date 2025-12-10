"""
Ring Component Analyzer
=======================

Identifies and classifies ring components using pattern matching
against known heterocycles and fused systems.

Reference: P-22, P-25
"""

from typing import List, Dict, Optional, Tuple, Set
from collections import Counter
import sys
import os

from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_structures import RingInfo


class ComponentMatcher:
    """
    Matches ring systems to known components using SMARTS patterns.
    
    Reference: P-22.2.1 (retained heterocycles), P-25.1 (fused hydrocarbons),
               P-25.2 (fused heterocycles)
    """
    
    # SMARTS patterns for component identification
    SMARTS_PATTERNS = {
        # =====================================================================
        # P-25.2.1 RETAINED FUSED HETEROCYCLES
        # =====================================================================
        
        # Tricyclic
        'phenazine': ['c1ccc2nc3ccccc3nc2c1'],
        'acridine': ['c1ccc2nc3ccccc3cc2c1'],
        'phenothiazine': ['c1ccc2Nc3ccccc3Sc2c1'],
        'phenoxazine': ['c1ccc2Nc3ccccc3Oc2c1'],
        'carbazole': ['c1ccc2[nH]c3ccccc3c2c1'],
        'xanthene': ['c1ccc2Oc3ccccc3Cc2c1'],
        
        # Bicyclic benzofused
        'quinoline': ['c1ccc2ncccc2c1'],
        'isoquinoline': ['c1ccc2ccncc2c1'],
        'quinazoline': ['c1ccc2ncncc2c1'],
        'quinoxaline': ['c1ccc2nccnc2c1'],
        'phthalazine': ['c1ccc2cnncc2c1'],
        'cinnoline': ['c1ccc2nncc2c1'],  # Note: different pattern
        'indole': ['c1ccc2[nH]ccc2c1', 'c1ccc2cc[nH]c2c1'],
        'benzofuran': ['c1ccc2occc2c1', 'c1ccc2ccoc2c1'],
        'benzothiophene': ['c1ccc2sccc2c1', 'c1ccc2ccsc2c1'],
        'benzimidazole': ['c1ccc2[nH]cnc2c1', 'c1ccc2nc[nH]c2c1'],
        'benzoxazole': ['c1ccc2ocnc2c1', 'c1ccc2ncoc2c1'],
        'benzothiazole': ['c1ccc2scnc2c1', 'c1ccc2ncsc2c1'],
        
        # Purine and pteridine
        'purine': ['c1ncc2[nH]cnc2n1', 'c1nc2[nH]cnc2c(=O)[nH]1'],
        'pteridine': ['c1cnc2nccnc2n1'],
        
        # =====================================================================
        # P-25.1.1 RETAINED FUSED HYDROCARBONS
        # =====================================================================
        'naphthalene': ['c1ccc2ccccc2c1'],
        'anthracene': ['c1ccc2cc3ccccc3cc2c1'],
        'phenanthrene': ['c1ccc2c(c1)ccc1ccccc12'],
        'pyrene': ['c1cc2ccc3cccc4ccc(c1)c2c34'],
        'azulene': ['c1ccc2cccc-2cc1'],  # 7-5 fused
        'fluorene': ['c1ccc2Cc3ccccc3c2c1'],
        'fluoranthene': ['c1ccc2c(c1)-c1cccc3cccc-2c13'],
        
        # =====================================================================
        # P-22.2.1 RETAINED MONOCYCLIC HETEROCYCLES
        # =====================================================================
        
        # 6-membered with 1 N
        'pyridine': ['n1ccccc1', 'c1ccncc1', 'c1cccnc1'],
        
        # 6-membered with 2 N
        'pyrazine': ['n1ccncc1', 'c1nccnc1'],
        'pyrimidine': ['n1cnccc1', 'c1ncncc1', 'c1ccncn1'],
        'pyridazine': ['n1ncccc1', 'c1ccnnc1'],
        
        # 6-membered with 3 N
        '1,3,5-triazine': ['n1cncnc1'],
        '1,2,4-triazine': ['n1ncncc1'],
        '1,2,3-triazine': ['n1nncc1'],
        
        # 5-membered with 1 heteroatom
        'furan': ['o1cccc1', 'c1ccoc1'],
        'thiophene': ['s1cccc1', 'c1ccsc1'],
        'pyrrole': ['[nH]1cccc1', 'c1cc[nH]c1'],
        'selenophene': ['[se]1cccc1'],
        'tellurophene': ['[te]1cccc1'],
        
        # 5-membered with 2 N
        'imidazole': ['n1ccnc1', 'c1nccn1', 'c1[nH]cnc1', 'c1nc[nH]c1'],
        'pyrazole': ['n1nccc1', 'c1ccnn1', 'c1cc[nH]n1'],
        
        # 5-membered with N and O
        'oxazole': ['o1ccnc1', 'n1ccoc1'],
        'isoxazole': ['o1nccc1', 'n1occc1'],
        
        # 5-membered with N and S
        'thiazole': ['s1ccnc1', 'n1ccsc1'],
        'isothiazole': ['s1nccc1', 'n1sccc1'],
        
        # 5-membered with 3 N
        '1,2,3-triazole': ['n1nncc1', 'c1nn[nH]c1'],
        '1,2,4-triazole': ['n1ncnc1', 'c1nc[nH]n1'],
        
        # 5-membered with 4 N
        'tetrazole': ['n1nnnc1', 'c1nnn[nH]1'],
        
        # =====================================================================
        # P-22.1 CARBOCYCLIC
        # =====================================================================
        'benzene': ['c1ccccc1'],
        'cyclopentadiene': ['C1=CC=CC1'],
        'cycloheptatriene': ['C1=CC=CC=CC1'],
    }
    
    # Fusion prefixes (P-25.3)
    FUSION_PREFIXES = {
        # Carbocyclic
        'benzene': 'benzo',
        'naphthalene': 'naphtho',
        'anthracene': 'anthro',
        'phenanthrene': 'phenanthro',
        'cyclopentadiene': 'cyclopenta',
        
        # Monocyclic heterocycles
        'pyridine': 'pyrido',
        'pyrazine': 'pyrazino',
        'pyrimidine': 'pyrimido',
        'pyridazine': 'pyridazino',
        'furan': 'furo',
        'thiophene': 'thieno',
        'pyrrole': 'pyrrolo',
        'imidazole': 'imidazo',
        'pyrazole': 'pyrazolo',
        'oxazole': 'oxazolo',
        'isoxazole': 'isoxazolo',
        'thiazole': 'thiazolo',
        'isothiazole': 'isothiazolo',
        '1,2,3-triazole': '[1,2,3]triazolo',
        '1,2,4-triazole': '[1,2,4]triazolo',
        'tetrazole': 'tetrazolo',
        '1,3,5-triazine': '[1,3,5]triazino',
        
        # Fused heterocycles
        'quinoline': 'quinolino',
        'isoquinoline': 'isoquinolino',
        'indole': 'indolo',
        'benzofuran': 'benzofuro',
        'benzothiophene': 'benzothieno',
        'benzimidazole': 'benzimidazolo',
        'purine': 'purino',
        'phenazine': 'phenazino',
        'acridine': 'acridino',
        'carbazole': 'carbazolo',
    }
    
    # Priority for base component selection (P-44.2, P-25.8)
    # Higher = more senior = preferred as base
    PRIORITIES = {
        # Tricyclic fused - highest
        'phenazine': 150,
        'acridine': 145,
        'phenothiazine': 144,
        'phenoxazine': 143,
        'carbazole': 142,
        'xanthene': 140,
        
        # Bicyclic fused
        'quinoline': 130,
        'isoquinoline': 130,
        'quinazoline': 128,
        'quinoxaline': 128,
        'phthalazine': 126,
        'cinnoline': 126,
        'indole': 125,
        'benzofuran': 124,
        'benzothiophene': 123,
        'benzimidazole': 127,
        'benzoxazole': 122,
        'benzothiazole': 121,
        'purine': 135,
        'pteridine': 134,
        
        # Fused hydrocarbons
        'naphthalene': 120,
        'anthracene': 119,
        'phenanthrene': 118,
        'azulene': 121,
        'pyrene': 117,
        'fluorene': 115,
        
        # 6-membered monocycles with N
        'pyridine': 100,
        'pyrazine': 99,
        'pyrimidine': 99,
        'pyridazine': 99,
        '1,3,5-triazine': 98,
        
        # 5-membered with 2+ N
        'imidazole': 90,
        'pyrazole': 89,
        '1,2,3-triazole': 88,
        '1,2,4-triazole': 88,
        'tetrazole': 87,
        
        # 5-membered with N + O/S
        'oxazole': 85,
        'isoxazole': 84,
        'thiazole': 83,
        'isothiazole': 82,
        
        # 5-membered with 1 heteroatom
        'furan': 70,
        'thiophene': 69,
        'pyrrole': 68,
        
        # Carbocyclic
        'benzene': 50,
        'cyclopentadiene': 40,
    }
    
    @classmethod
    def identify_component(cls, mol: Chem.Mol, ring_atoms: Set[int]) -> Optional[str]:
        """
        Identify a ring component by matching against known patterns.
        
        Args:
            mol: RDKit molecule
            ring_atoms: Set of atom indices in the ring/system
        
        Returns:
            Component name or None
        """
        for name, patterns in cls.SMARTS_PATTERNS.items():
            for smarts in patterns:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is None:
                    continue
                
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    match_set = set(match)
                    # Check if this match corresponds to our ring
                    overlap = len(match_set & ring_atoms)
                    if overlap >= len(ring_atoms) * 0.7:  # 70% overlap
                        return name
        
        return None
    
    @classmethod
    def find_all_components(cls, mol: Chem.Mol, 
                           rings: List[RingInfo]) -> List[Dict]:
        """
        Find all identifiable components in a molecule.
        
        Returns list of dicts with:
        - name: component name
        - atoms: set of atom indices
        - priority: seniority value
        - fusion_prefix: prefix for fusion names
        """
        components = []
        used_atoms = set()
        
        # Get all ring atoms
        all_ring_atoms = set()
        for ring in rings:
            all_ring_atoms.update(ring.atoms)
        
        # First pass: large fused systems
        large_systems = [
            'phenazine', 'acridine', 'phenothiazine', 'phenoxazine',
            'carbazole', 'xanthene', 'quinoline', 'isoquinoline',
            'quinazoline', 'quinoxaline', 'indole', 'benzofuran',
            'benzothiophene', 'benzimidazole', 'naphthalene',
            'anthracene', 'phenanthrene', 'purine', 'pteridine'
        ]
        
        for name in large_systems:
            patterns = cls.SMARTS_PATTERNS.get(name, [])
            for smarts in patterns:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is None:
                    continue
                
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    match_set = set(match)
                    
                    # Must be part of ring system
                    if not match_set.issubset(all_ring_atoms):
                        continue
                    
                    # Check overlap with already used atoms
                    overlap = len(match_set & used_atoms)
                    if overlap < len(match_set) * 0.3:
                        components.append({
                            'name': name,
                            'atoms': match_set,
                            'priority': cls.PRIORITIES.get(name, 0),
                            'fusion_prefix': cls.FUSION_PREFIXES.get(name, name + 'o')
                        })
                        used_atoms.update(match_set)
                        break
        
        # Second pass: monocyclic rings
        for ring in rings:
            ring_set = set(ring.atoms)
            
            # Skip if mostly used
            if len(ring_set & used_atoms) >= len(ring_set) * 0.7:
                continue
            
            # Try to identify
            for name, patterns in cls.SMARTS_PATTERNS.items():
                if name in large_systems:
                    continue
                
                for smarts in patterns:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern is None:
                        continue
                    
                    matches = mol.GetSubstructMatches(pattern)
                    for match in matches:
                        match_set = set(match)
                        if match_set & ring_set:
                            # Found a match including this ring
                            if len(match_set & used_atoms) < len(match_set) * 0.5:
                                components.append({
                                    'name': name,
                                    'atoms': match_set,
                                    'priority': cls.PRIORITIES.get(name, 0),
                                    'fusion_prefix': cls.FUSION_PREFIXES.get(name, name + 'o')
                                })
                                used_atoms.update(match_set)
                                break
                    else:
                        continue
                    break
        
        return components
    
    @classmethod
    def get_base_component(cls, components: List[Dict]) -> Optional[Dict]:
        """
        Select the base (most senior) component.
        
        Reference: P-25.3.4.2.1, P-44.2
        """
        if not components:
            return None
        
        # Sort by priority (highest first), then by size (largest first)
        sorted_components = sorted(
            components,
            key=lambda c: (-c['priority'], -len(c['atoms']))
        )
        
        return sorted_components[0]
