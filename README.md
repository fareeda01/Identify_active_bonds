# Filtering Significant Bond Length Changes in a Reaction Step (Active Region)
This Python script analyzes structural changes during a chemical reaction by comparing 
interatomic bond lengths across three molecular geometries:
- Reactant (REA)
- Transition State (TS)
- Product (PRO)
  
## Features

Parses atomic coordinates from .xyz files Identifies bonded atom pairs using covalent radii and a customizable distance toleranceComputes and compares bond lengths across REA, TS, and PRO structuresDetects significant bond length changes (e.g., bond breaking/forming)Outputs results to:
bond_length_changes.csv — all bondssignificant_bond_changes.csv — filtered by a significance threshold which can be changed depending on the problem


### Dependencies:  
numpy  Standard Python libraries (csv, itertools)
