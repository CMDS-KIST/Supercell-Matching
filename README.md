# supercell_matching

Find common supercells of given two lattices where ab planes of two lattices are similar (e.g. square-square, hexagonal-hexagonal)

## Requirements
pymatgen

numpy

sympy

## Usage
1. Create a pre-calculated table to find common supercells with a maximum index and the constant D determined by the symmetry (e.g. square: D = -1, hexagonal: D = -3)
```
python Zquad_maketable.py --max_index {number} D
```
2. Find common supercells of two lattice cif files from the pre-calculated table
```
python find_common_supercell_main.py --tol_strain {float} --output_type {cif|latex} --target_angles_file {none|filename} --max_index {number} --tol_angle {float} Top_cif_filename Bottom_cif_filename 
```
