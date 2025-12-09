"""
IUPAC Namer - Main Interface
============================

Main class for generating IUPAC names from molecular structures.

Usage:
    from rdkit import Chem
    from iupac_nomenclature import IUPACNamer

    mol = Chem.MolFromSmiles('CCO')
    namer = IUPACNamer()
    name = namer.name(mol)
    print(name)  # ethanol
"""

from typing import Optional, List, Dict, Union
from rdkit import Chem

from iupac_naming.analyzers import MolecularAnalyzer, ComponentMatcher
from iupac_naming.builders import NameBuilder, FusedNameBuilder
from iupac_naming.data_structures import MoleculeAnalysis
from iupac_naming.rules import (
    RETAINED_HETEROCYCLES, RETAINED_FUSED_HETEROCYCLES,
    RETAINED_FUSED_HYDROCARBONS, RETAINED_FUNCTIONAL_PARENTS
)


class IUPACNamer:
    """
    Main class for IUPAC nomenclature generation.

    This class provides the primary interface for converting molecular
    structures to IUPAC names following the 2013 Blue Book recommendations.

    Attributes:
        verbose: Print debug information during naming
        prefer_retained: Use retained (trivial) names when available
        prefer_pin: Generate Preferred IUPAC Names (PINs)

    Example:
        >>> namer = IUPACNamer(verbose=True)
        >>> mol = Chem.MolFromSmiles('c1ccccc1')
        >>> namer.name(mol)
        'benzene'
    """

    def __init__(self,
                 verbose: bool = False,
                 prefer_retained: bool = True,
                 prefer_pin: bool = True):
        """
        Initialize the IUPAC namer.

        Args:
            verbose: Print debug information
            prefer_retained: Use retained names when available
            prefer_pin: Generate Preferred IUPAC Names
        """
        self.verbose = verbose
        self.prefer_retained = prefer_retained
        self.prefer_pin = prefer_pin

    def _log(self, msg: str):
        """Print debug message if verbose."""
        if self.verbose:
            print(f"[IUPACNamer] {msg}")

    # =========================================================================
    # MAIN NAMING METHODS
    # =========================================================================

    def name(self, mol: Chem.Mol) -> str:
        """
        Generate IUPAC name for a molecule.

        This is the main entry point for nomenclature generation.

        Args:
            mol: RDKit Mol object

        Returns:
            IUPAC name string

        Raises:
            ValueError: If molecule is invalid
        """
        if mol is None:
            raise ValueError("Invalid molecule (None)")

        self._log(f"Naming molecule with {mol.GetNumAtoms()} atoms")

        # Step 1: Quick check for simple retained names
        if self.prefer_retained:
            retained = self._check_retained_name(mol)
            if retained:
                self._log(f"Using retained name: {retained}")
                return retained

        # Step 2: Analyze molecular structure
        analyzer = MolecularAnalyzer(mol)
        analysis = analyzer.analyze()

        self._log(f"Structure type: {analysis.parent_type}")
        self._log(f"Rings: {len(analysis.rings)}")
        self._log(f"Functional groups: {len(analysis.characteristic_groups)}")

        # Step 3: Route to appropriate naming strategy
        if analysis.parent_type == 'fused':
            return self._name_fused_system(mol, analysis)
        elif analysis.parent_type == 'spiro':
            return self._name_spiro_system(mol, analysis)
        elif analysis.parent_type == 'monocyclic':
            return self._name_monocyclic(mol, analysis)
        elif analysis.parent_type == 'acyclic':
            return self._name_acyclic(mol, analysis)
        else:
            return self._name_generic(mol, analysis)

    def name_from_smiles(self, smiles: str) -> str:
        """
        Generate IUPAC name from SMILES string.

        Args:
            smiles: SMILES string

        Returns:
            IUPAC name
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return self.name(mol)

    def name_from_molblock(self, molblock: str) -> str:
        """
        Generate IUPAC name from MOL block.

        Args:
            molblock: MOL format string (V2000 or V3000)

        Returns:
            IUPAC name
        """
        mol = Chem.MolFromMolBlock(molblock)
        if mol is None:
            raise ValueError("Invalid MOL block")
        return self.name(mol)

    # =========================================================================
    # RETAINED NAME CHECKING
    # =========================================================================

    def _check_retained_name(self, mol: Chem.Mol) -> Optional[str]:
        """
        Check if molecule matches a retained name.

        Reference: P-12.3, P-22.2.1, P-25.1.1, P-25.2.1
        """
        # Check simple heterocycles
        for name, info in RETAINED_HETEROCYCLES.items():
            pattern_smarts = self._get_smarts_for_retained(name)
            if pattern_smarts:
                pattern = Chem.MolFromSmarts(pattern_smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    if mol.GetNumAtoms() == pattern.GetNumAtoms():
                        return name

        # Check fused hydrocarbons
        for name in RETAINED_FUSED_HYDROCARBONS:
            patterns = ComponentMatcher.SMARTS_PATTERNS.get(name, [])
            for smarts in patterns:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    if mol.GetNumAtoms() == pattern.GetNumAtoms():
                        return name

        # Check fused heterocycles
        for name in RETAINED_FUSED_HETEROCYCLES:
            patterns = ComponentMatcher.SMARTS_PATTERNS.get(name, [])
            for smarts in patterns:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    if mol.GetNumAtoms() == pattern.GetNumAtoms():
                        return name

        return None

    def _get_smarts_for_retained(self, name: str) -> Optional[str]:
        """Get SMARTS pattern for a retained name."""
        patterns = ComponentMatcher.SMARTS_PATTERNS.get(name, [])
        return patterns[0] if patterns else None

    # =========================================================================
    # SPECIFIC NAMING METHODS
    # =========================================================================

    def _name_fused_system(self, mol: Chem.Mol,
                          analysis: MoleculeAnalysis) -> str:
        """
        Name a fused ring system.

        Reference: P-25
        """
        self._log("Naming fused ring system")

        # Find components
        components = ComponentMatcher.find_all_components(mol, analysis.rings)

        self._log(f"Found {len(components)} components")
        for comp in components:
            self._log(f"  - {comp['name']} (priority {comp['priority']})")

        # Use specialized fused name builder with mol object
        builder = FusedNameBuilder(analysis, components, mol=mol, verbose=self.verbose)
        return builder.build_name()

    def _name_spiro_system(self, mol: Chem.Mol,
                          analysis: MoleculeAnalysis) -> str:
        """
        Name a spiro compound.

        Reference: P-24
        """
        self._log("Naming spiro system")
        builder = NameBuilder(analysis, self.verbose)
        return builder.build_name()

    def _name_monocyclic(self, mol: Chem.Mol,
                        analysis: MoleculeAnalysis) -> str:
        """
        Name a monocyclic compound.

        Reference: P-22
        """
        self._log("Naming monocyclic compound")

        # Try to find retained name for the ring
        if analysis.rings:
            ring = analysis.rings[0]
            ring_atoms = set(ring.atoms)

            for name, patterns in ComponentMatcher.SMARTS_PATTERNS.items():
                for smarts in patterns:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern:
                        matches = mol.GetSubstructMatches(pattern)
                        for match in matches:
                            if set(match) == ring_atoms:
                                self._log(f"Ring identified as: {name}")

                                # Check for substituents
                                if mol.GetNumAtoms() == len(ring.atoms):
                                    return name

                                # Has substituents - build full name
                                builder = NameBuilder(analysis, self.verbose)
                                return builder.build_name()

        builder = NameBuilder(analysis, self.verbose)
        return builder.build_name()

    def _name_acyclic(self, mol: Chem.Mol,
                     analysis: MoleculeAnalysis) -> str:
        """
        Name an acyclic compound.

        Reference: P-21, P-59.1.1
        """
        self._log("Naming acyclic compound")
        builder = NameBuilder(analysis, self.verbose)
        return builder.build_name()

    def _name_generic(self, mol: Chem.Mol,
                     analysis: MoleculeAnalysis) -> str:
        """Name using generic approach."""
        self._log("Using generic naming")
        builder = NameBuilder(analysis, self.verbose)
        return builder.build_name()

    # =========================================================================
    # UTILITY METHODS
    # =========================================================================

    def analyze(self, mol: Chem.Mol) -> MoleculeAnalysis:
        """
        Get structural analysis without naming.

        Useful for debugging or understanding molecular structure.

        Args:
            mol: RDKit Mol object

        Returns:
            MoleculeAnalysis object
        """
        analyzer = MolecularAnalyzer(mol)
        return analyzer.analyze()

    def get_components(self, mol: Chem.Mol) -> List[Dict]:
        """
        Get identified ring components.

        Args:
            mol: RDKit Mol object

        Returns:
            List of component dictionaries
        """
        analyzer = MolecularAnalyzer(mol)
        rings = analyzer.get_rings()
        return ComponentMatcher.find_all_components(mol, rings)