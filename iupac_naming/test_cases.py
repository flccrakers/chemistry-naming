"""
IUPAC Nomenclature Test Cases
=============================

Test cases for validating IUPAC nomenclature generation.
All tests are based on IUPAC 2013 Blue Book recommendations.
"""
from iupac_naming.namer import IUPACNamer
from iupac_naming.tests.main_test import test_cases_smiles, test_case_mol_V3000
from iupac_naming.tests.enovalys_test_part1 import enovalys_test_part1

test_cases_smiles = test_cases_smiles

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
