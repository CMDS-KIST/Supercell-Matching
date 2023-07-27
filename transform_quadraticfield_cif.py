from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.lattice import Lattice
import os
import argparse
import numpy as np

default_workdir = 'transform_quadfield'

def str_transformed(dsquare: int):
    return f'transform_quadfield_{dsquare}'
def get_transformed_quadfield_surface_for_heterostructure(initial_structure, transform_operator, dsquare, is_top = True):
    transformed_structure = initial_structure.copy()
    if transform_operator.shape == (2,):
        m1 = transform_operator[0]
        m2 = transform_operator[1]
        transformed_structure.make_supercell([[m1, m2, 0], [dsquare*m2, m1, 0], [0, 0, 1]])
    elif transform_operator.shape == (2, 2):
        m11 = transform_operator[0][0]
        m12 = transform_operator[0][1]
        m21 = transform_operator[1][0]
        m22 = transform_operator[1][1]
        transformed_structure.make_supercell([[m11, m12, 0], [m21, m22, 0], [0, 0, 1]])
    
    a, b, c = transformed_structure.lattice.matrix
    transformed_vacuum_lattice = [a, b, c*2]
    transformed_vacuum_fcoord = transformed_structure.frac_coords.copy()
    transformed_vacuum_fcoord[:, 2] = (transformed_vacuum_fcoord[:, 2] - np.floor(transformed_vacuum_fcoord[:, 2] + 1e-5))/2.0
    if not is_top:
        transformed_vacuum_fcoord[:, 2] = transformed_vacuum_fcoord[:, 2] + 0.5
    transformed_vacuum_added_structure = Structure(transformed_vacuum_lattice, transformed_structure.species_and_occu, transformed_vacuum_fcoord)
    return transformed_vacuum_added_structure, transformed_structure

def export_transformed_quadfield_surface_for_heterostructure(initial_structure_filename, transform_operator, dsquare, is_top = True, workdir = default_workdir):
    if not os.path.isdir(workdir):
        os.mkdir(workdir)

    initial_cif_parser = CifParser(initial_structure_filename)
    initial_structure = initial_cif_parser.get_structures(primitive=False)[0]
    transformed_structure, transformed_structure2 = get_transformed_quadfield_surface_for_heterostructure(initial_structure, transform_operator, dsquare, is_top)
    initial_basename = os.path.basename(initial_structure_filename)
    initial_filename_without_ext = os.path.splitext(initial_basename)[0]
    if transform_operator.shape == (2,):
        m1 = transform_operator[0]
        m2 = transform_operator[1]        
        transformed_filename = os.path.join(default_workdir, f'{initial_filename_without_ext}_{m1}I_{m2}sqrt{dsquare}.cif')
        transformed_filename2 = os.path.join(default_workdir, f'{initial_filename_without_ext}_{m1}I_{m2}sqrt{dsquare}_novacuum.cif')
    elif transform_operator.shape == (2, 2):
        m11 = transform_operator[0][0]
        m12 = transform_operator[0][1]
        m21 = transform_operator[1][0]
        m22 = transform_operator[1][1]
        transformed_filename = os.path.join(default_workdir, f'{initial_filename_without_ext}_{m11}_{m12}_{m21}_{m22}.cif')
        transformed_filename2 = os.path.join(default_workdir, f'{initial_filename_without_ext}_{m11}_{m12}_{m21}_{m22}_novacuum.cif')
    transformed_vacuum_added_cif_writer = CifWriter(transformed_structure)
    transformed_vacuum_added_cif_writer.write_file(transformed_filename)
    transformed_cif_writer = CifWriter(transformed_structure2)
    transformed_cif_writer.write_file(transformed_filename2)
    return transformed_structure

def export_heterostructure(top_structure_filename, transform_operator1, bottom_structure_filename, transform_operator2, dsquare, workdir = default_workdir):
    if not os.path.isdir(workdir):
        os.mkdir(workdir)
 
    top_transformed = export_transformed_quadfield_surface_for_heterostructure(top_structure_filename, transform_operator1, dsquare, workdir = workdir)
    bottom_transformed = export_transformed_quadfield_surface_for_heterostructure(bottom_structure_filename, transform_operator2, dsquare, is_top=False, workdir = workdir)

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
    if transform_operator1.shape == (2,) and transform_operator2.shape == (2,):
        m1 = transform_operator1[0]
        m2 = transform_operator1[1]
        n1 = transform_operator2[0]
        n2 = transform_operator2[1]
        heterostructure_filename = os.path.join(workdir, f'{top_filename_without_ext}_{m1}I+{m2}sqrt{dsquare}_{bottom_filename_without_ext}_{n1}I+{n2}sqrt{dsquare}.cif')
        heterostructure_reduced_filename = os.path.join(workdir, f'{top_filename_without_ext}_{m1}I+{m2}sqrt{dsquare}_{bottom_filename_without_ext}_{n1}I+{n2}sqrt{dsquare}_reduced.cif')
    else:
        heterostructure_filename = os.path.join(workdir, f'{top_filename_without_ext}_T1_{bottom_filename_without_ext}_T2.cif')
        heterostructure_reduced_filename = os.path.join(workdir, f'{top_filename_without_ext}_T1_{bottom_filename_without_ext}_T2_reduced.cif')
    heterostructure_cif_writer.write_file(heterostructure_filename)

    heterostructure_reduced = heterostructure.get_reduced_structure(reduction_algo='LLL')
    
    heterostructure_reduced_cif_writer = CifWriter(heterostructure_reduced)
    heterostructure_reduced_cif_writer.write_file(heterostructure_reduced_filename)

    return heterostructure, heterostructure_filename

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Top_cif_filename', type=str, help='The cif file name for top surface.')
    parser.add_argument('m1', type=int, help='m1 for m1*I + m2*sqrt(d) quadratic field transformation for top surface')
    parser.add_argument('m2', type=int, help='m2 for m1*I + m2*sqrt(d) quadratic field transformation for top surface')
    parser.add_argument('Bottom_cif_filename', type=str, help='The cif file name for bottom surface.')
    parser.add_argument('n1', type=int, help='n1 for n1*I + n2*sqrt(d) quadratic field transformation for top surface')
    parser.add_argument('n2', type=int, help='n2 for n1*I + n2*sqrt(d) quadratic field transformation for top surface')
    parser.add_argument('dsquare', type=int, help='d square spanning quadratic field representing lattice')

    args = parser.parse_args()

    print(f'{args.Top_cif_filename} is transformed by {args.m1}*I + {args.m2}*Omega')
    print(f'{args.Bottom_cif_filename} is transformed by {args.n1}*I + {args.n2}*Omega')
    heterostructure, filename = export_heterostructure(args.Top_cif_filename, np.array([args.m1, args.m2]), args.Bottom_cif_filename, np.array([args.n1, args.n2]), args.dsquare)
    sa_prec1 = SpacegroupAnalyzer(heterostructure, symprec=0.1, angle_tolerance=10)
    heterostructure_prim = sa_prec1.find_primitive()
    print(f'a: {heterostructure_prim.lattice.a}, b: {heterostructure_prim.lattice.b}, c: {heterostructure_prim.lattice.c}, alpha: {heterostructure_prim.lattice.alpha}, beta: {heterostructure_prim.lattice.beta}, gamma: {heterostructure_prim.lattice.gamma}')

    sa_prec2 = SpacegroupAnalyzer(heterostructure, symprec=0.01)
    heterostructure_prim = sa_prec2.find_primitive()
    print(f'a: {heterostructure_prim.lattice.a}, b: {heterostructure_prim.lattice.b}, c: {heterostructure_prim.lattice.c}, alpha: {heterostructure_prim.lattice.alpha}, beta: {heterostructure_prim.lattice.beta}, gamma: {heterostructure_prim.lattice.gamma}')

    sa_prec3 = SpacegroupAnalyzer(heterostructure, symprec=0.001)
    heterostructure_prim = sa_prec3.find_primitive()
    print(f'a: {heterostructure_prim.lattice.a}, b: {heterostructure_prim.lattice.b}, c: {heterostructure_prim.lattice.c}, alpha: {heterostructure_prim.lattice.alpha}, beta: {heterostructure_prim.lattice.beta}, gamma: {heterostructure_prim.lattice.gamma}')
    CifWriter(heterostructure_prim).write_file(f'{filename}.prim.cif')

    #smaller_supercell = export_heterostructure(args.Top_cif_filename, np.array([[1, 1], [-2, 1]]), args.Bottom_cif_filename, np.array([[2, 0], [1, 1]]), args.dsquare)
    #smaller_supercell.make_supercell([[1, 0, 0], [-1, 2, 0], [0, 0, 1]])
    #print(f'a:{smaller_supercell.lattice.a}, b:{smaller_supercell.lattice.b}, gamma:{smaller_supercell.lattice.gamma}')
    #reduced_smaller_supercell = smaller_supercell.get_reduced_structure(reduction_algo='LLL')
    #print(f'a:{reduced_smaller_supercell.lattice.a}, b:{reduced_smaller_supercell.lattice.b}, gamma:{reduced_smaller_supercell.lattice.gamma}')

if __name__ == "__main__":
    main()