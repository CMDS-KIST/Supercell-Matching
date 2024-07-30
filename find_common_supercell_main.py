import argparse
import enum
import os
import cmath, math
import numpy as np
from print_supercell_text import GaussInt_square, EisensteinInt_hex
import transform_hexagonal_cif
import transform_quadraticfield_cif
import bisect

from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

work_dir = 'heterostructure_temp'
float_equal_tol = 1e-4

def determinant_lattice_shape(l: Lattice):
    return l.b/l.a*(math.cos(l.gamma/180*math.pi) + 1j*math.sin(l.gamma/180*math.pi))

def find_ratio(target_value, ratio_table, d:int, tol=1e-2, target_angles=[]):
    if d > 0:
        print('Quadratic integer is real - it is not 2D lattice!')
        exit(1)

    def make_omega_matrix(m11, m12, d:int):
        if d != -3:
            m11i_plus_m12omega = np.array([[m11, m12], [d*m12, m11]])
        elif d == -3:
            m11i_plus_m12omega = np.array([[m11, m12], [-m12, m11-m12]])
        return m11i_plus_m12omega

    def make_return_dict(i):
        m = ratio_table[i][0][0] + ratio_table[i][0][1]*omega
        n = ratio_table[i][0][2] + ratio_table[i][0][3]*omega
        angle = cmath.phase(m/n)/math.pi*180

        return_dict = {'Transform_a': make_omega_matrix(ratio_table[i][0][0], ratio_table[i][0][1], d),
            'Transform_b': make_omega_matrix(ratio_table[i][0][2], ratio_table[i][0][3], d),
            'Quad_int_a': [ratio_table[i][0][0], ratio_table[i][0][1]],
            'Quad_int_b': [ratio_table[i][0][2], ratio_table[i][0][3]],
            'angle': angle,
            'table_ratio': ratio_table[i][1]}
        return return_dict

    if d != -3:
        omega = cmath.sqrt(d)
    elif d == -3:
        omega = (-1 + cmath.sqrt(d))/2.0

    sorted_ratio_value = [value[1] for value in ratio_table]
    n_left = bisect.bisect_left(sorted_ratio_value, target_value)

    i = n_left
    return_dict_list = []
    return_dict_list.append(make_return_dict(i))
    i = i + 1
    while (abs(sorted_ratio_value[n_left] - sorted_ratio_value[i]) < tol*target_value):
        if not target_angles:
            return_dict_list.append(make_return_dict(i))
        else:
            pass
        i = i + 1
    
    return return_dict_list

def find_ratio_by_quad_pair(target_value, ratio_table, d:int, tol=1e-2, target_angles = []):
    if d > 0:
        print('Quadratic integer is real - it is not 2D lattice!')
        exit(1)

    if d != -3:
        omega = cmath.sqrt(d)
    elif d == -3:
        omega = (-1 + cmath.sqrt(d))/2.0

    sorted_ratio_value = [value[1] for value in ratio_table]
    n_left = bisect.bisect_left(sorted_ratio_value, target_value)

    i = n_left
    return_dict_list = []
    if d != -3 and d != -1:
        return_dict_list.append([ratio_table[i][0][0] + ratio_table[i][0][1]*omega, ratio_table[i][0][2] + ratio_table[i][0][3]*omega, 0])
    elif d == -3:
        return_dict_list.append([EisensteinInt_hex(ratio_table[i][0][0], ratio_table[i][0][1]), EisensteinInt_hex(ratio_table[i][0][2], ratio_table[i][0][3]), 0])
    elif d == -1:
        return_dict_list.append([GaussInt_square(ratio_table[i][0][0], ratio_table[i][0][1]), GaussInt_square(ratio_table[i][0][2], ratio_table[i][0][3]), 0])

    i = i + 1
    i_neg = i - 1
    previous_index = 0
    candadite_list = []
    for j in range(len(target_angles)):
        candadite_list.append([])
    while (abs(sorted_ratio_value[n_left] - sorted_ratio_value[i]) < tol*target_value):
        if d != -3 and d != -1:
            layer1_quadint = ratio_table[i][0][0] + ratio_table[i][0][1]*omega
            layer2_quadint = ratio_table[i][0][2] + ratio_table[i][0][3]*omega
            layer1_angle = cmath.phase(layer1_quadint)
            layer2_angle = cmath.phase(layer2_quadint)
        elif d == -3:
            layer1_quadint = EisensteinInt_hex(ratio_table[i][0][0], ratio_table[i][0][1])
            layer2_quadint = EisensteinInt_hex(ratio_table[i][0][2], ratio_table[i][0][3])
            layer1_angle = layer1_quadint.arg()
            layer2_angle = layer2_quadint.arg()
        elif d == -1:
            layer1_quadint = GaussInt_square(ratio_table[i][0][0], ratio_table[i][0][1])
            layer2_quadint = GaussInt_square(ratio_table[i][0][2], ratio_table[i][0][3])
            layer1_angle = layer1_quadint.arg()
            layer2_angle = layer2_quadint.arg()

        rotation_angle = (layer1_angle - layer2_angle)/math.pi*180
        strain = (sorted_ratio_value[i] - target_value)/target_value
        solution = [layer1_quadint, layer2_quadint, strain]
        #print(rotation_angle)
        if not target_angles:
            return_dict_list.append(solution)
        else:
            index_for_target_angle = np.where(np.isclose(abs(rotation_angle), target_angles, rtol=0.001, atol=0.01))[0]
            if index_for_target_angle.size != 0:
                found_index = index_for_target_angle[0]
                candadite_list[found_index].append(solution)
        i = i + 1

    while (abs(sorted_ratio_value[n_left] - sorted_ratio_value[i_neg]) < tol*target_value):
        if d != -3 and d != -1:
            layer1_quadint = ratio_table[i_neg][0][0] + ratio_table[i_neg][0][1]*omega
            layer2_quadint = ratio_table[i_neg][0][2] + ratio_table[i_neg][0][3]*omega
            layer1_angle = cmath.phase(layer1_quadint)
            layer2_angle = cmath.phase(layer2_quadint)
        elif d == -3:
            layer1_quadint = EisensteinInt_hex(ratio_table[i_neg][0][0], ratio_table[i_neg][0][1])
            layer2_quadint = EisensteinInt_hex(ratio_table[i][0][2], ratio_table[i_neg][0][3])
            layer1_angle = layer1_quadint.arg()
            layer2_angle = layer2_quadint.arg()
        elif d == -1:
            layer1_quadint = GaussInt_square(ratio_table[i_neg][0][0], ratio_table[i_neg][0][1])
            layer2_quadint = GaussInt_square(ratio_table[i_neg][0][2], ratio_table[i_neg][0][3])
            layer1_angle = layer1_quadint.arg()
            layer2_angle = layer2_quadint.arg()

        rotation_angle = (layer1_angle - layer2_angle)/math.pi*180
        strain = (sorted_ratio_value[i_neg] - target_value)/target_value
        solution = [layer1_quadint, layer2_quadint, strain]
        #print(rotation_angle)
        if not target_angles:
            return_dict_list.append(solution)
        else:
            index_for_target_angle = np.where(np.isclose(abs(rotation_angle), target_angles, rtol=0.001, atol=0.01))[0]
            if index_for_target_angle.size != 0:
                found_index = index_for_target_angle[0]
                candadite_list[found_index].append(solution)
        i_neg = i_neg - 1
    
    if candadite_list:
        for candadites in candadite_list:
            if candadites:
                if d != -3 and d != -1:
                    min_index = min(enumerate(map(lambda item:abs(item[0]) + abs(item[1]), candadites)), key=lambda x:x[1])[0]
                else:
                    min_index = min(enumerate(map(lambda item:item[0].norm() + item[1].norm(), candadites)), key=lambda x:x[1])[0]
                return_dict_list.append(candadites[min_index])

    if d == -3 or d == -1:
        #return_dict_list.sort(key=lambda x: abs(cmath.phase(x[0].complex_form()) - cmath.phase(x[1].complex_form())))
        return_dict_list.sort(key=lambda x: abs(x[0].arg() - x[1].arg()))
    else:
        return_dict_list.sort(key=lambda x: abs(cmath.phase(x[0]) - cmath.phase(x[1])))
    return return_dict_list

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

    if abs(d_top**2 + d_top + 1) < float_equal_tol:
        import eisenstein_list_pos_arg as eisenstein_list
        D = -3
        ratio_table = eisenstein_list.sorted_dict
    elif d_top.real < 1e-4:
        D = -1*round(d_top.imag**2)
        if abs(d_top.imag**2 - 1) < float_equal_tol:
            import gaussian_list_pos_arg
            ratio_table = gaussian_list_pos_arg.sorted_dict
        elif abs(d_top.imag**2 - 2) < float_equal_tol:
            import Zsqrt2i_list
            ratio_table = Zsqrt2i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 5) < float_equal_tol:
            import Zsqrt5i_list 
            ratio_table = Zsqrt5i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 6) < float_equal_tol:
            import Zsqrt6i_list
            ratio_table = Zsqrt6i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 7) < float_equal_tol:
            import Zsqrt7i_list
            ratio_table = Zsqrt7i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 10) < float_equal_tol:
            import Zsqrt10i_list
            ratio_table = Zsqrt10i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 11) < float_equal_tol:
            import Zsqrt11i_list
            ratio_table = Zsqrt11i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 13) < float_equal_tol:
            import Zsqrt13i_list
            ratio_table = Zsqrt13i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 14) < float_equal_tol:
            import Zsqrt14i_list
            ratio_table = Zsqrt14i_list.sorted_ratio_dict
        elif abs(d_top.imag**2 - 15) < float_equal_tol:
            import Zsqrt15i_list
            ratio_table = Zsqrt15i_list.sorted_ratio_dict
        else:
            print('This rectangular lattice is not supported')
            exit(1)
    else:
        print('This oblique lattice is not supported.')
        exit(1)
    
    target_ratio = Bottom_structure.lattice.a/Top_structure.lattice.a
    Top_cif_filename = args.Top_cif_filename
    Bottom_cif_filename = args.Bottom_cif_filename
    if target_ratio < 1:
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
        roots = find_ratio(target_ratio, ratio_table, D, args.tol_strain, target_angle_list)
        for root in roots:
            if D == -3:
                heterostructure, filename = transform_hexagonal_cif.export_heterostructure(Top_cif_filename, root['Quad_int_a'][0], root['Quad_int_a'][1], Bottom_cif_filename, root['Quad_int_b'][0], root['Quad_int_b'][1], workdir = work_dir)
            else:
                heterostructure, filename = transform_quadraticfield_cif.export_heterostructure(Top_cif_filename, np.array(root['Quad_int_a']), Bottom_cif_filename, np.array(root['Quad_int_b']), D, workdir = work_dir)

            sa_prec1 = SpacegroupAnalyzer(heterostructure, symprec=0.1, angle_tolerance=10)
            heterostructure_prim = sa_prec1.find_primitive()

            base_filename = os.path.basename(filename)
            base_splited_filename = os.path.splitext(base_filename)[0]
            output_filename = f'{base_splited_filename}_prim.cif'
            CifWriter(heterostructure_prim).write_file(output_filename)

            print_common_supercell_transform_dict(root, target_ratio)
    elif args.output_type == 'latex':
        import print_supercell_text
        roots = find_ratio_by_quad_pair(target_ratio, ratio_table, D, args.tol_strain, target_angle_list)
        for root in roots:
            print(print_supercell_text.common_supercell_wood_notations(root[0:2], root[2], D))

if __name__ == '__main__':
    import time
    btime = time.time()
    main()
    etime = time.time()
    print(etime - btime)