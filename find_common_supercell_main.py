import argparse
import enum
import os
import cmath, math
import numpy as np
from print_supercell_text import GaussInt_square, EisensteinInt_hex
import transform_hexagonal_cif
import transform_quadraticfield_cif
import bisect
import pickle
from quadraticinteger import DedekindDomainInt

from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

work_dir = 'heterostructure_temp'
float_equal_tol = 1e-4

def determinant_lattice_shape(l: Lattice):
    return l.b/l.a*(math.cos(l.gamma/180*math.pi) + 1j*math.sin(l.gamma/180*math.pi))

def find_ratio(target_value, ratio_table, quad_d:int, tol=1e-2, target_angles=[], tol_angle = 0.01):
    if quad_d > 0:
        print('This quadratic integer is real - it is not 2D lattice!')
        exit(1)

    def make_omega_matrix(m11:int, m12:int, quad_d:int):
        if quad_d % 4 != 1:
            m11i_plus_m12omega = np.array([[m11, m12], [quad_d*m12, m11]])
        else:
            m11i_plus_m12omega = np.array([[m11, m12], [(quad_d - 1)//4*m12, m11 - m12]])
        return m11i_plus_m12omega

    def make_return_dict(i):
        return_dict = {'Transform_a': make_omega_matrix(ratio_table[i][0][0], ratio_table[i][0][1], quad_d),
            'Transform_b': make_omega_matrix(ratio_table[i][0][2], ratio_table[i][0][3], quad_d),
            'Quad_int_a': [ratio_table[i][0][0], ratio_table[i][0][1]],
            'Quad_int_b': [ratio_table[i][0][2], ratio_table[i][0][3]],
            'angle': ratio_table[i][1][1],
            'table_ratio': ratio_table[i][1][0]}
        return return_dict

    sorted_ratio_value = [value[1][0] for value in ratio_table]
    target_lower_bound = target_value*(1.0 - tol)
    target_upper_bound = target_value*(1.0 + tol)
    n_left = bisect.bisect_left(sorted_ratio_value, target_lower_bound)
    n_right = bisect.bisect_right(sorted_ratio_value, target_upper_bound)

    r_filtered_dict_list = [make_return_dict(i) for i in range(n_left, n_right)]
    if not target_angles:
        return_dict_list = r_filtered_dict_list
    else:
        angle_tol = tol_angle
        angle_filtered_dict_list = [d for d in r_filtered_dict_list if any(abs(d['angle'] - target_angle) < angle_tol for target_angle in target_angles)]
        if quad_d == -1:
            angle_filtered_dict_list += [d for d in r_filtered_dict_list if any(abs(d['angle'] - (90 - target_angle)) < angle_tol for target_angle in target_angles)]
        elif quad_d == -3:
            angle_filtered_dict_list += [d for d in r_filtered_dict_list if any(abs(d['angle'] - (60 - target_angle)) < angle_tol for target_angle in target_angles)]
        return_dict_list = angle_filtered_dict_list
    
    return return_dict_list

def find_ratio_by_quad_pair(target_value, ratio_table, quad_d:int, tol=1e-2, target_angles = [], tol_angle=0.01):
    if quad_d > 0:
        print('Quadratic integer is real - it is not 2D lattice!')
        exit(1)

    sorted_ratio_value = [value[1][0] for value in ratio_table]
    target_lower_bound = target_value*(1.0 - tol)
    target_upper_bound = target_value*(1.0 + tol)
    n_left = bisect.bisect_left(sorted_ratio_value, target_lower_bound)
    n_right = bisect.bisect_right(sorted_ratio_value, target_upper_bound)

    r_filtered_quad_pair = [[DedekindDomainInt(ratio_table[i][0][0], ratio_table[i][0][1], quad_d), DedekindDomainInt(ratio_table[i][0][2], ratio_table[i][0][3], quad_d), ratio_table[i][1][1], (ratio_table[i][1][0] - target_value)/target_value] for i in range(n_left, n_right)]
    if not target_angles:
        return_quad_pair = r_filtered_quad_pair
    else:
        angle_tol = tol_angle
        angle_filtered_quad_pair = [pair for pair in r_filtered_quad_pair if any(abs(pair[2] - target_angle) < angle_tol for target_angle in target_angles)]
        return_quad_pair = angle_filtered_quad_pair
        if quad_d == -1:
            angle_filtered_quad_pair += [pair for pair in r_filtered_quad_pair if any(abs(pair[2] - (90 - target_angle)) < angle_tol for target_angle in target_angles)]
        elif quad_d == -3:
            angle_filtered_quad_pair += [pair for pair in r_filtered_quad_pair if any(abs(pair[2] - (60 - target_angle)) < angle_tol for target_angle in target_angles)]

        return_quad_pair.sort(key=lambda x: x[2])
    return return_quad_pair

def print_common_supercell_transform_dict(transform_dict, test_r):
    transform_dict_str = f'[{transform_dict["Transform_a"][0][0]: 03d} {transform_dict["Transform_a"][0][1]: 03d}]l_a = exp(i*{transform_dict["angle"]:0.02f}/180*pi)[{transform_dict["Transform_b"][0][0]: 03d} {transform_dict["Transform_b"][0][1]: 03d}]({transform_dict["table_ratio"]/test_r})l_b ({transform_dict["table_ratio"]:.4f})\n[{transform_dict["Transform_a"][1][0]: 03d} {transform_dict["Transform_a"][1][1]: 03d}]                          [{transform_dict["Transform_b"][0][0]: 03d} {transform_dict["Transform_b"][1][1]: 03d}]'
    print(transform_dict_str) 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Top_cif_filename', type=str, help='The cif file name for top surface.')
    parser.add_argument('Bottom_cif_filename', type=str, help='The cif file name for bottom surface.')
    parser.add_argument('--tol_strain', type=float, default=0.01, help='Tolerance to strain (default is 0.01)')
    parser.add_argument('--output_type', type=str, default='cif', help='Type of output (cif, latex)')
    parser.add_argument('--target_angles_file', type=str, default='', help='File name describing rotation angles of the heterostructures to finding.')
    parser.add_argument('--max_index', type=int, default=60, help='Maximum limit of the coefficients of the quadratic integers')
    parser.add_argument('--tol_angle', type=float, default=0.01, help='Tolerance to angle (default is 0.01)')
    
    args = parser.parse_args()
    if not os.path.isfile(args.Top_cif_filename):
        print('No top cif file')
        exit(1)
    if not os.path.isfile(args.Bottom_cif_filename):
        print('No bottom cif file')
        exit(1)

    Top_cif_parser = CifParser(args.Top_cif_filename)
    Bottom_cif_parser = CifParser(args.Bottom_cif_filename)
    Top_structure = Top_cif_parser.get_structures(primitive=False)[0]
    Bottom_structure = Bottom_cif_parser.get_structures(primitive=False)[0]

    d_top = determinant_lattice_shape(Top_structure.lattice)
    d_bottom = determinant_lattice_shape(Bottom_structure.lattice)

    if abs(d_top - d_bottom) > float_equal_tol:
        print('Shapes of top and bottom are different.')
        exit(1)

    if abs(d_top.real + 0.5) < float_equal_tol:
        quad_d = -round(4*d_top.imag**2)
    elif abs(d_top.real) < float_equal_tol:
        quad_d = -round(d_top.imag**2)
    else:
        print('This lattice is not supported.')
        exit(1)
    
    filename_pickle = f'OZsqrt{-quad_d}i_list_{args.max_index}.pickle'
    if os.path.isfile(filename_pickle):
        with open(filename_pickle, 'rb') as f:
            ratio_table = pickle.load(f)
    else:
        print('There is no table to find common supercells. Please make a table by executing Zquad_maketable.py.')
        exit(1)
    
    target_ratio = Bottom_structure.lattice.a/Top_structure.lattice.a
    Top_cif_filename = args.Top_cif_filename
    Bottom_cif_filename = args.Bottom_cif_filename
    if target_ratio < 1: # If the bottom cell is smaller than the top cell, we swap the top and the bottom cells.
        temp_structure = Top_structure.copy()
        Top_structure = Bottom_structure.copy()
        Bottom_structure = temp_structure.copy()
        target_ratio = Bottom_structure.lattice.a/Top_structure.lattice.a
        Top_cif_filename = args.Bottom_cif_filename
        Bottom_cif_filename = args.Top_cif_filename

    target_angle_list = []
    if args.target_angles_file != '':
        if os.path.isfile(args.target_angles_file):
            with open(args.target_angles_file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    target_angle_list.append(float(line))
        else:
            print('Angle list file does not exist.')
            exit(1)

    if args.output_type == 'cif':
        roots = find_ratio(target_ratio, ratio_table, quad_d, args.tol_strain, target_angle_list, args.tol_angle)
        for root in roots:
            heterostructure, filename = transform_quadraticfield_cif.export_heterostructure(Top_cif_filename, np.array(root['Quad_int_a']), Bottom_cif_filename, np.array(root['Quad_int_b']), quad_d, workdir = work_dir)

            #sa_prec1 = SpacegroupAnalyzer(heterostructure, symprec=0.1, angle_tolerance=10)
            #heterostructure_prim = sa_prec1.find_primitive()

            #base_filename = os.path.basename(filename)
            #base_splited_filename = os.path.splitext(base_filename)[0]
            #output_filename = f'{base_splited_filename}_prim.cif'
            #CifWriter(heterostructure_prim).write_file(output_filename)

            print_common_supercell_transform_dict(root, target_ratio)
    elif args.output_type == 'latex':
        import print_supercell_text
        roots = find_ratio_by_quad_pair(target_ratio, ratio_table, quad_d, args.tol_strain, target_angle_list, args.tol_angle)
        for root in roots:
            print(print_supercell_text.common_supercell_wood_notations(root[0:3], root[3], quad_d))

if __name__ == '__main__':
    import time
    btime = time.time()
    main()
    etime = time.time()
    print(etime - btime)