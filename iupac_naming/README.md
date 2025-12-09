# chemistry-naming
# IUPAC Nomenclature System

Implementation of IUPAC 2013 Blue Book nomenclature rules for systematic naming of organic compounds.

## Architecture

```
iupac_nomenclature/
├── __init__.py              # Main package exports
├── namer.py                 # IUPACNamer - main entry point
├── data_structures.py       # Data classes for molecular representation
├── test_cases.py            # Test cases (SMILES and V3000 MOL)
│
├── analyzers/               # Molecular structure analysis
│   ├── __init__.py
│   ├── molecular_analyzer.py    # MolecularAnalyzer - extracts structure info
│   └── component_matcher.py     # ComponentMatcher - identifies ring components
│
├── builders/                # Name construction
│   ├── __init__.py
│   ├── name_builder.py          # NameBuilder, FusedNameBuilder
│   └── (fusion_builder.py)      # TODO: Specialized fusion name builder
│
├── rules/                   # IUPAC Blue Book rules implementation
│   ├── __init__.py
│   ├── p1_general.py            # P-1: General principles, locants, prefixes
│   ├── p2_parent_hydrides.py    # P-2: Parent hydrides, rings, fused systems
│   ├── p3_characteristic_groups.py  # P-3: Functional groups, suffixes
│   ├── p4_name_construction.py  # P-4: Name construction rules, seniority
│   ├── (p5_selection.py)        # TODO: P-5 Selecting PINs
│   ├── (p6_applications.py)     # TODO: P-6 Specific compound classes
│   ├── (p7_radicals.py)         # TODO: P-7 Radicals, ions
│   ├── (p8_isotopes.py)         # TODO: P-8 Isotopic modifications
│   └── (p9_stereochemistry.py)  # TODO: P-9 Configuration/conformation
│
├── utils/                   # Helper functions
│   └── __init__.py              # Locant formatting, elision, etc.
│
└── data/                    # Reference data (future)
    └── (element_data.py)        # TODO: Element properties, priorities
```

## Blue Book Section Mapping

Each function/constant is traceable to specific Blue Book sections:

| Section | Topic | Implementation |
|---------|-------|----------------|
| P-12 | Name types | `rules/p1_general.py::NameType` |
| P-13 | Operations | `rules/p1_general.py::NomenclatureOperation` |
| P-14.1 | Bonding numbers | `rules/p1_general.py::get_standard_bonding_number()` |
| P-14.2 | Multiplicative prefixes | `rules/p1_general.py::get_multiplicative_prefix()` |
| P-14.3 | Locants | `rules/p1_general.py::format_locants()` |
| P-14.5 | Alphanumerical order | `rules/p1_general.py::alphanumerical_sort_key()` |
| P-14.7 | Indicated H | `rules/p1_general.py::format_indicated_hydrogen()` |
| P-16 | Name writing | `rules/p1_general.py::NameWriter` |
| P-21 | Acyclic hydrides | `rules/p2_parent_hydrides.py::get_acyclic_hydride_name()` |
| P-22.1 | Monocyclic | `rules/p2_parent_hydrides.py::get_cycloalkane_name()` |
| P-22.2 | Hantzsch-Widman | `rules/p2_parent_hydrides.py::get_hantzsch_widman_name()` |
| P-22.2.1 | Retained heterocycles | `rules/p2_parent_hydrides.py::RETAINED_HETEROCYCLES` |
| P-25.1 | Fused hydrocarbons | `rules/p2_parent_hydrides.py::RETAINED_FUSED_HYDROCARBONS` |
| P-25.2 | Fused heterocycles | `rules/p2_parent_hydrides.py::RETAINED_FUSED_HETEROCYCLES` |
| P-25.3 | Fusion prefixes | `rules/p2_parent_hydrides.py::FUSION_PREFIXES` |
| P-25.4 | Bridges | `rules/p2_parent_hydrides.py::BIVALENT_BRIDGES` |
| P-28 | Ring assemblies | `rules/p2_parent_hydrides.py::get_ring_assembly_name()` |
| P-31 | Unsaturation | `rules/p3_characteristic_groups.py::format_unsaturation()` |
| P-33 | Functional suffixes | `rules/p3_characteristic_groups.py::FUNCTIONAL_SUFFIXES` |
| P-34 | Retained parents | `rules/p3_characteristic_groups.py::RETAINED_FUNCTIONAL_PARENTS` |
| P-35 | Substituent prefixes | `rules/p3_characteristic_groups.py::CHARACTERISTIC_GROUP_PREFIXES` |
| P-41 | Class seniority | `rules/p4_name_construction.py::CLASS_SENIORITY` |
| P-44.2 | Ring seniority | `rules/p4_name_construction.py::compare_ring_seniority()` |
| P-44.3 | Chain selection | `rules/p4_name_construction.py::select_principal_chain()` |
| P-44.4 | Low locant rule | `rules/p4_name_construction.py::apply_low_locant_rule()` |

## Usage

```python
from rdkit import Chem
from iupac_nomenclature import IUPACNamer

namer = IUPACNamer(verbose=True)

# From SMILES
name = namer.name_from_smiles('c1ccccc1')  # -> benzene

# From RDKit Mol
mol = Chem.MolFromSmiles('c1ccc2nc3ccccc3nc2c1')
name = namer.name(mol)  # -> phenazine

# Get structural analysis
analysis = namer.analyze(mol)
print(f"Parent type: {analysis.parent_type}")
print(f"Rings: {len(analysis.rings)}")
```

## Current Status

### Working:
- ✅ Retained names (pyridine, benzene, naphthalene, phenazine, etc.)
- ✅ Monocyclic heterocycles (P-22.2.1)
- ✅ Fused hydrocarbons (P-25.1.1)
- ✅ Fused heterocycles (P-25.2.1)
- ✅ Component identification in fused systems
- ✅ Epoxy bridge detection
- ✅ Basic ring analysis

### In Progress:
- 🔄 Fusion locant calculation (P-25.3.2, P-25.3.3)
- 🔄 Peripheral numbering system
- 🔄 Acyclic compound naming
- 🔄 Substituent handling

### TODO:
- ⬜ Complete fusion locant system with primes
- ⬜ Spiro nomenclature (P-24)
- ⬜ von Baeyer bridged systems (P-23)
- ⬜ Ring assemblies (P-28)
- ⬜ Stereochemistry (P-9)
- ⬜ Isotopes (P-8)
- ⬜ Radicals and ions (P-7)

## Test Results

```
Extended test suite: 18/19 passed (95%)

Retained names working:
- pyridine, pyrazine, pyrimidine, pyridazine
- furan, thiophene, pyrrole
- imidazole, oxazole, thiazole
- naphthalene, anthracene
- quinoline, isoquinoline, indole
- benzofuran, benzothiophene
- benzimidazole, quinoxaline, quinazoline
- phenazine
```

## Complex Example

Target: `9,12-epoxypyrido[1'',2'':1',2']imidazo[4',5':5,6]pyrazino[2,3-b]phenazine`

Current output: `?,?-epoxypyridoimidazopyrazinophenazine`

Analysis:
- ✅ Components correctly identified: phenazine, pyrazine, imidazole, pyridine
- ✅ Epoxy bridge detected at correct positions
- ✅ Fusion prefixes correct (pyrido, imidazo, pyrazino)
- ⬜ Missing: Fusion locants with primes
- ⬜ Missing: Bridge locants

## References

- IUPAC Recommendations 2013 (Blue Book)
- Part 1: P-1 to P-4 (General principles, Parent hydrides, Groups, Construction)
- Part 2: P-5 to P-7 (Selection, Applications, Radicals)  
- Part 3: P-8 to P-9 (Isotopes, Stereochemistry)
