from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import os
import argparse
import numpy as np

default_work_dir = "transform_hexagonal"

def get_transformed_hexagonal_structure(initial_structure : Structure, m1, m2):
    transformed_structure = initial_structure.copy()
    transformed_structure.make_supercell([[m1, m2, 0], [-m2, m1 - m2, 0], [0, 0, 1]])
    for site in transformed_structure.sites:
        site.frac_coords[2] = site.frac_coords[2] - np.floor(site.frac_coords[2])

    return transformed_structure

def match_planar_coordinates_both(top_structure, bottom_structure):
    matched_top_structure = top_structure.copy()
    matched_bottom_structure = bottom_structure.copy()

    top_a = matched_top_structure.lattice.a
    bottom_a = matched_bottom_structure.lattice.a
    average_a = (top_a + bottom_a)/2.0
    top_average_ratio = average_a/top_a
    bottom_average_ratio = average_a/bottom_a

    top_strain_value = top_average_ratio - 1.0
    matched_top_structure.apply_strain([top_strain_value, top_strain_value, 0])
    bottom_strain_value = bottom_average_ratio - 1.0
    matched_bottom_structure.apply_strain([bottom_strain_value, bottom_strain_value, 0])

    return matched_top_structure, matched_bottom_structure

def match_planar_coordinates_top_to_bottom(top_structure, bottom_structure):
    matched_top_structure = top_structure.copy()
    matched_bottom_structure = bottom_structure.copy()

    top_a = matched_top_structure.lattice.a
    bottom_a = matched_bottom_structure.lattice.a
    top_bottom_ratio = bottom_a/top_a

    top_strain_value = top_bottom_ratio - 1.0
    matched_top_structure.apply_strain([top_strain_value, top_strain_value, 0])

    return matched_top_structure

def get_out_of_plain_coords_in_heterostructure(top_structure, bottom_structure, vacuum = 15, distance = 0.0):
    c_top = top_structure.lattice.c
    c_bottom = bottom_structure.lattice.c
    total_c = c_top + c_bottom + vacuum + distance

    new_top_lattice_a, new_top_lattice_b, new_top_lattice_c = top_structure.lattice.matrix
    new_top_lattice_c = new_top_lattice_c/c_top*total_c
    new_bottom_lattice_a, new_bottom_lattice_b, new_bottom_lattice_c = bottom_structure.lattice.matrix
    new_bottom_lattice_c = new_bottom_lattice_c/c_bottom*total_c

    new_top_frac_coord = top_structure.frac_coords
    new_top_frac_coord[:, 2] = new_top_frac_coord[:, 2]*c_top/total_c

    new_top_structure = Structure([new_top_lattice_a, new_top_lattice_b, new_top_lattice_c], top_structure.species, new_top_frac_coord)
    for site in new_top_structure:
        site.c = site.c - np.floor(site.c) + c_bottom/total_c + distance/total_c
    new_bottom_structure = Structure([new_bottom_lattice_a, new_bottom_lattice_b, new_bottom_lattice_c], bottom_structure.species, bottom_structure.cart_coords, coords_are_cartesian=True)

    return new_top_structure, new_bottom_structure

def concat_structures(top_structure, bottom_structure):
    concatonated_structure = bottom_structure.copy()
    for site in top_structure:
        concatonated_structure.append(site.species, site.frac_coords)
    return concatonated_structure

def get_transformed_hexagonal_surface_for_heterostructure(initial_structure, m1, m2, is_top = True):
    transformed_structure = initial_structure.copy()
    transformed_structure.make_supercell([[m1, m2, 0], [-m2, m1 - m2, 0], [0, 0, 1]])
    
    a, b, c = transformed_structure.lattice.matrix
    transformed_vacuum_lattice = [a, b, c*2]
    transformed_vacuum_fcoord = transformed_structure.frac_coords.copy()
    transformed_vacuum_fcoord[:, 2] = (transformed_vacuum_fcoord[:, 2] - np.floor(transformed_vacuum_fcoord[:, 2] + 1e-5))/2.0
    if not is_top:
        transformed_vacuum_fcoord[:, 2] = transformed_vacuum_fcoord[:, 2] + 0.5
    transformed_vacuum_added_structure = Structure(transformed_vacuum_lattice, transformed_structure.species_and_occu, transformed_vacuum_fcoord)
    return transformed_vacuum_added_structure

def export_structure(basename, export_structure, output_dir = '', export_type = 'cif'):
    if 'cif' in export_type:   
        exporter = CifWriter(export_structure)
        filename = os.path.join(output_dir, f'{basename}.cif')
        exporter.write_file(filename)
    if 'vasp' in export_type:
        exporter = Poscar(export_structure)
        filename = os.path.join(output_dir, f'{basename}.vasp')
        exporter.write_file(filename)
    if 'vdir' in export_type:
        exporter = Poscar(export_structure)
        if not os.path.isdir(basename):
            os.mkdir(basename)
        filename = os.path.join(output_dir, basename, 'POSCAR')
        exporter.write_file(filename)

def export_transformed_hexagonal_surface_for_heterostructure(initial_structure_filename, m1, m2, is_top = True, work_dir = default_work_dir, export_type = 'cif'):
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)

    initial_cif_parser = CifParser(initial_structure_filename)
    initial_structure = initial_cif_parser.get_structures(primitive=False)[0]

    transformed_structure = get_transformed_hexagonal_surface_for_heterostructure(initial_structure, m1, m2, is_top)
    initial_basename = os.path.basename(initial_structure_filename)
    initial_filename_without_ext = os.path.splitext(initial_basename)[0]
    if 'cif' in export_type:   
        transformed_vacuum_added_cif_writer = CifWriter(transformed_structure)
        transformed_filename = os.path.join(work_dir, f'{initial_filename_without_ext}_{m1}I_{m2}Omega.cif')
        transformed_vacuum_added_cif_writer.write_file(transformed_filename)
    if 'vasp' in export_type:
        transformed_vacuum_added_writer = Poscar(transformed_structure)
        transformed_filename = os.path.join(work_dir, f'{initial_filename_without_ext}_{m1}I_{m2}Omega.poscar.vasp')
        transformed_vacuum_added_writer.write_file(transformed_filename)
    return transformed_structure

def export_heterostructure(top_structure_filename, m1, m2, bottom_structure_filename, n1, n2, work_dir = default_work_dir):
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)

    top_transformed = export_transformed_hexagonal_surface_for_heterostructure(top_structure_filename, m1, m2, work_dir=work_dir)
    bottom_transformed = export_transformed_hexagonal_surface_for_heterostructure(bottom_structure_filename, n1, n2, is_top=False, work_dir=work_dir)

    top_a = top_transformed.lattice.a
    bottom_a = bottom_transformed.lattice.a
    average_a = (top_a + bottom_a)/2.0
    top_average_ratio = average_a/top_a
    bottom_average_ratio = average_a/bottom_a

    top_strain_value = top_average_ratio - 1.0
    top_transformed.apply_strain([top_strain_value, top_strain_value, 0])
    bottom_strain_value = bottom_average_ratio - 1.0
    bottom_transformed.apply_strain([bottom_strain_value, bottom_strain_value, 0])
    
    heterostructure = Structure(top_transformed.lattice, top_transformed.species_and_occu + bottom_transformed.species_and_occu, np.concatenate((top_transformed.frac_coords, bottom_transformed.frac_coords), axis=0))
    heterostructure_cif_writer = CifWriter(heterostructure)

    top_basename = os.path.basename(top_structure_filename)
    top_filename_without_ext = os.path.splitext(top_basename)[0]
    
    bottom_basename = os.path.basename(bottom_structure_filename)
    bottom_filename_without_ext = os.path.splitext(bottom_basename)[0]
    heterostructure_filename = os.path.join(work_dir, f'{top_filename_without_ext}_{m1}I+{m2}Omega_{bottom_filename_without_ext}_{n1}I+{n2}Omega.cif')
    heterostructure_cif_writer.write_file(heterostructure_filename)

    return heterostructure, heterostructure_filename

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Top_cif_filename', type=str, help='The cif file name for top surface.')
    parser.add_argument('m1', type=int, help='m1 for m1*I + m2*Omega hexagonal transformation for top surface')
    parser.add_argument('m2', type=int, help='m2 for m1*I + m2*Omega hexagonal transformation for top surface')
    parser.add_argument('Bottom_cif_filename', type=str, help='The cif file name for bottom surface.')
    parser.add_argument('n1', type=int, help='n1 for n1*I + n2*Omega hexagonal transformation for top surface')
    parser.add_argument('n2', type=int, help='n2 for n1*I + n2*Omega hexagonal transformation for top surface')

    args = parser.parse_args()

    print(f'{args.Top_cif_filename} is transformed by {args.m1}*I + {args.m2}*Omega')
    print(f'{args.Bottom_cif_filename} is transformed by {args.n1}*I + {args.n2}*Omega')
    export_heterostructure(args.Top_cif_filename, args.m1, args.m2, args.Bottom_cif_filename, args.n1, args.n2)

if __name__ == "__main__":
    main()