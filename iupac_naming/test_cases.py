"""
IUPAC Nomenclature Test Cases
=============================

Test cases for validating IUPAC nomenclature generation.
All tests are based on IUPAC 2013 Blue Book recommendations.
"""
from iupac_naming import IUPACNamer

# =============================================================================
# SIMPLE SMILES TEST CASES
# =============================================================================

test_cases_smiles = [
    # P-22.2.1 Retained monocyclic heterocycles
    ("n1ccccc1", "pyridine"),
    ("n1ccncc1", "pyrazine"),
    ("c1c[nH]cn1", "imidazole"),
    ("o1cccc1", "furan"),
    ("s1cccc1", "thiophene"),
    ("[nH]1cccc1", "pyrrole"),
    ("n1cnccc1", "pyrimidine"),
    ("n1ncccc1", "pyridazine"),
    ("c1cc[nH]n1", "pyrazole"),
    ("o1ccnc1", "oxazole"),
    ("s1ccnc1", "thiazole"),
    
    # P-25.1.1 Retained fused hydrocarbons
    ("c1ccc2ccccc2c1", "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1", "anthracene"),
    ("c1ccc2c(c1)ccc1ccccc12", "phenanthrene"),
    
    # P-25.2.1 Retained fused heterocycles
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccc2ccncc2c1", "isoquinoline"),
    ("c1ccc2[nH]ccc2c1", "indole"),
    ("c1ccc2occc2c1", "benzofuran"),
    ("c1ccc2sccc2c1", "benzothiophene"),
    ("c1ccc2nc3ccccc3nc2c1", "phenazine"),
    ("c1ccc2[nH]cnc2c1", "benzimidazole"),
    ("c1ccc2nccnc2c1", "quinoxaline"),
    ("c1ccc2ncncc2c1", "quinazoline"),
    
    # P-22.1 Carbocycles
    ("c1ccccc1", "benzene"),
    ("C1CCCCC1", "cyclohexane"),
    ("C1CCCC1", "cyclopentane"),
    ("C1CCC1", "cyclobutane"),
    ("C1CC1", "cyclopropane"),
    
    # P-21.2.1 Acyclic hydrocarbons
    ("C", "methane"),
    ("CC", "ethane"),
    ("CCC", "propane"),
    ("CCCC", "butane"),
    ("CCCCC", "pentane"),
    ("C=C", "ethene"),
    ("C#C", "ethyne"),
    ("CC=CC", "but-2-ene"),
    
    # P-63 Alcohols
    ("CO", "methanol"),
    ("CCO", "ethanol"),
    ("CCCO", "propan-1-ol"),
    ("CC(O)C", "propan-2-ol"),
    
    # P-64 Ketones
    ("CC(=O)C", "acetone"),
    ("CCC(=O)C", "butan-2-one"),
    
    # P-65 Carboxylic acids
    ("C(=O)O", "formic acid"),
    ("CC(=O)O", "acetic acid"),
    ("CCC(=O)O", "propanoic acid"),
    
    # P-66 Aldehydes
    ("C=O", "formaldehyde"),
    ("CC=O", "acetaldehyde"),
    
    # P-62 Amines
    ("CN", "methylamine"),
    ("CCN", "ethylamine"),
    ("c1ccccc1N", "aniline"),
    ## EDGE CASES AND COMPLEX EXAMPLES
    # Single atoms
    ("C", "methane"),
    ("N", "azane"),  # ammonia
    ("O", "oxidane"),  # water

    # Simple chains
    ("CCCCCCCCCC", "decane"),
    ("CCCCCCCCCCCCCCCCCCCCC", "henicosane"),  # 21 carbons

    # Highly substituted
    ("CC(C)(C)C", "2,2-dimethylpropane"),  # neopentane
    ("CC(C)C(C)C", "2,3-dimethylbutane"),

    # Multiple functional groups
    ("OCC(O)CO", "propane-1,2,3-triol"),  # glycerol
    ("NC(=O)c1ccccc1", "benzamide"),
    # Substituted heterocycles
    ("Cc1ncccc1", "2-methylpyridine"),  # 2-picoline
    ("Cc1ccccn1", "2-methylpyridine"),  # alternate SMILES
    ("c1ccc(O)cc1", "phenol"),
    ("c1ccc(N)cc1", "aniline"),
    ("c1ccc(C(=O)O)cc1", "benzoic acid"),

    # Multiple substituents
    ("Cc1ccc(O)cc1", "4-methylphenol"),  # p-cresol
    ("Oc1ccc(O)cc1", "benzene-1,4-diol"),  # hydroquinone
    ("Nc1ccc(N)cc1", "benzene-1,4-diamine"),

    # Fused with substituents
    ("Oc1ccc2ccccc2c1", "naphthalen-2-ol"),  # 2-naphthol
    ("Cc1ccc2ccccc2c1", "2-methylnaphthalene"),

    # Polycyclic
    ("c1cc2ccc3cccc4ccc(c1)c2c34", "pyrene"),
    ("c1ccc2c(c1)ccc1c2ccc2ccccc21", "chrysene"),

]

test_case_mol_V3000 = [
    # ---------------------------------------------------------------------
    # Complex fused heterocyclic with epoxy bridge
    # Reference: P-25 (fused systems), P-25.4 (bridges)
    # ---------------------------------------------------------------------
    ("""9,12-epoxypyrido[1'',2'':1',2']imidazo[4',5':5,6]pyrazino[2,3-b]phenazine
  Mrv2120 12092515372D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 26 32 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.2263 0.9287 0 0
M  V30 2 C 0.3927 -0.6023 0 0
M  V30 3 C 1.8018 -1.2237 0 0
M  V30 4 C 3.0444 -0.3141 0 0
M  V30 5 N 4.4535 -0.9354 0 0
M  V30 6 C 5.6962 -0.0258 0 0
M  V30 7 C 7.1052 -0.6472 0 0
M  V30 8 C 6.7724 2.4148 0 0
M  V30 9 C 5.5297 1.5051 0 0
M  V30 10 N 4.1207 2.1265 0 0
M  V30 11 C 2.878 1.2169 0 0
M  V30 12 C 1.4689 1.8383 0 0
M  V30 13 C 8.1815 1.7934 0 0
M  V30 14 C 8.3479 0.2624 0 0
M  V30 15 N 9.757 -0.359 0 0
M  V30 16 N 9.4241 2.703 0 0
M  V30 17 N 12.2378 2.713 0 0
M  V30 18 C 10.8332 2.0816 0 0
M  V30 19 C 10.9996 0.5507 0 0
M  V30 20 N 12.5071 0.2358 0 0
M  V30 21 C 13.2724 1.5722 0 0
M  V30 22 C 14.8124 1.5777 0 0
M  V30 23 C 15.5871 0.2468 0 0
M  V30 24 C 14.8218 -1.0896 0 0
M  V30 25 C 13.2818 -1.0951 0 0
M  V30 26 O 1.2526 0.2657 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 14
M  V30 8 2 13 8
M  V30 9 1 8 9
M  V30 10 1 6 9
M  V30 11 2 9 10
M  V30 12 1 10 11
M  V30 13 1 4 11
M  V30 14 2 11 12
M  V30 15 1 1 12
M  V30 16 1 13 14
M  V30 17 1 14 15
M  V30 18 2 15 19
M  V30 19 2 18 16
M  V30 20 1 13 16
M  V30 21 2 21 17
M  V30 22 1 17 18
M  V30 23 1 18 19
M  V30 24 1 20 19
M  V30 25 1 20 21
M  V30 26 1 21 22
M  V30 27 2 22 23
M  V30 28 1 23 24
M  V30 29 2 24 25
M  V30 30 1 20 25
M  V30 31 1 12 26
M  V30 32 1 3 26
M  V30 END BOND
M  V30 END CTAB
M  END
""", "9,12-epoxypyrido[1'',2'':1',2']imidazo[4',5':5,6]pyrazino[2,3-b]phenazine"),

    # ---------------------------------------------------------------------
    # Trispiro compound
    # Reference: P-24 (spiro systems)
    # ---------------------------------------------------------------------
    ("""Trispiro compound
  Mrv2120 12092516512D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 23 27 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 7.3487 1.102 0 0
M  V30 2 C 8.88 0.9383 0 0
M  V30 3 C 9.5039 -0.4697 0 0
M  V30 4 C 8.5965 -1.714 0 0
M  V30 5 C 7.0652 -1.5503 0 0
M  V30 6 C 6.4413 -0.1423 0 0
M  V30 7 C 5.4089 -1.285 0 0
M  V30 8 C 5.6737 1.1927 0 0
M  V30 9 C 4.1668 0.8751 0 0
M  V30 10 C 4.0031 -0.6561 0 0
M  V30 11 C 1.3508 -0.3727 0 0 CFG=1
M  V30 12 C 2.5951 -1.2801 0 0
M  V30 13 O 2.9225 1.7825 0 0
M  V30 14 C 1.5145 1.1586 0 0
M  V30 15 C -0.5038 -2.0368 0 0
M  V30 16 C 1.028 -1.8785 0 0
M  V30 17 C 0.0185 0.3996 0 0
M  V30 18 C -1.1277 -0.6288 0 0
M  V30 19 C -3.666 -1.4486 0 0
M  V30 20 C -2.1602 -1.7715 0 0
M  V30 21 C -1.6011 0.8366 0 0
M  V30 22 C -3.1068 1.1595 0 0
M  V30 23 C -4.1393 0.0168 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 1 6
M  V30 7 2 8 9
M  V30 8 1 9 10
M  V30 9 2 7 10
M  V30 10 1 7 6
M  V30 11 1 6 8
M  V30 12 1 11 12
M  V30 13 1 13 14
M  V30 14 1 11 14
M  V30 15 1 12 10
M  V30 16 1 9 13
M  V30 17 1 15 16
M  V30 18 1 17 18
M  V30 19 1 15 18
M  V30 20 1 11 16 CFG=3
M  V30 21 1 11 17
M  V30 22 1 19 20
M  V30 23 1 21 22
M  V30 24 1 22 23
M  V30 25 1 19 23
M  V30 26 1 20 18
M  V30 27 1 18 21
M  V30 END BOND
M  V30 END CTAB
M  END
""", "2''H,4''H-trispiro[cyclohexane-1,1'-cyclopentane-3',3''-cyclopenta[b]pyran-6'',1'''-cyclohexane]"),
]


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_tests():
    """Run all test cases and report results."""
    from rdkit import Chem

    
    namer = IUPACNamer(verbose=False)
    
    print("=" * 70)
    print("IUPAC NOMENCLATURE TEST SUITE")
    print("=" * 70)
    
    passed = 0
    failed = 0
    
    for smiles, expected in test_cases_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  INVALID SMILES: {smiles}")
            failed += 1
            continue
        
        try:
            result = namer.name(mol)
            # Accept equivalent names
            acceptable = [expected.lower()]
            if expected.lower() == 'acetone':
                acceptable.append('propan-2-one')
            
            if result.lower() in acceptable:
                print(f"  PASS: {smiles:25} -> {result}")
                passed += 1
            else:
                print(f"  FAIL: {smiles:25} -> {result} (expected: {expected})")
                failed += 1
        except Exception as e:
            print(f"  ERROR: {smiles:25} -> {e}")
            failed += 1

    for molblock, expected in test_case_mol_V3000:
        mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=True)
        if mol is None:
            print(f"  INVALID MOLBLOCK")
            failed += 1
            continue

        try:
            result = namer.name(mol)
            if result.lower() == expected.lower():
                print(f"  PASS: V3000 MOLBLOCK -> {result}")
                passed += 1
            else:
                print(f"  FAIL: V3000 MOLBLOCK -> {result} (expected: {expected})")
                failed += 1
        except Exception as e:
            print(f"  ERROR: V3000 MOLBLOCK -> {e}")
            failed += 1
    
    print("=" * 70)
    print(f"Results: {passed} passed, {failed} failed")
    print(f"Success rate: {100*passed/(passed+failed):.1f}%")
    

if __name__ == "__main__":
    run_tests()
