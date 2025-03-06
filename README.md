# Supercell-Matching

Find common supercells of given two lattices where ab planes of two lattices are similar (e.g. square-square, hexagonal-hexagonal)

## Requirements
* Python
* Pymatgen
* Numpy
* Sympy

## Usage
1. Create a pre-calculated table to find common supercells with a maximum index and the constant D determined by the symmetry (e.g. square: D = -1, hexagonal: D = -3)
```
python Zquad_maketable.py --max_index {number} D
```
2. Find common supercells of two lattice cif files from the pre-calculated table
```
python find_common_supercell_main.py --tol_strain {float} --output_type {cif|latex} --target_angles_file {none|filename} --max_index {number} --tol_angle {float} Top_cif_filename Bottom_cif_filename 
```

## Citation
If you use this code, please cite our paper:
```
@article{https://doi.org/10.1002/smtd.202400579,
author = {Lee, Weon-Gyu and Lee, Jung-Hoon},
title = {A Deterministic Method to Construct a Common Supercell Between Two Similar Crystalline Surfaces},
journal = {Small Methods},
volume = {8},
number = {12},
pages = {2400579},
doi = {https://doi.org/10.1002/smtd.202400579},
year = {2024}
}
```
