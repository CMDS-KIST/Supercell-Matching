import math, cmath, bisect
import operator
import eisenstein
import numpy as np
import os
import sys

def make_eisenstein_supercell_dict(max_index):
    omega = -0.5+cmath.sqrt(-3)/2
    #omega = cmath.sqrt(-1)
    ratio_dict = {}
    for i in range(0, max_index):
        for j in range(i+1, max_index):
            for k in range(0, max_index):
                for l in range(k+1, max_index):
                    if i == 0 and k == 0:
                        m_eisenstein = eisenstein.EisensteinInt(j, i)
                        n_eisenstein = eisenstein.EisensteinInt(l, k)
                        m_n_gcd = m_eisenstein.gcd(n_eisenstein)
                        _, m_arg = m_eisenstein.polar_form()
                        _, n_arg = n_eisenstein.polar_form()
                        arg_diff = m_arg - n_arg
                        if m_n_gcd.is_unit() and (i*i - i*j + j*j > k*k - k*l + l*l or (i*i - i*j + j*j == k*k - k*l + l*l and arg_diff >= 0)):
                            ratio_dict[(j, i, l, k)] = math.sqrt((i*i - i*j + j*j)/(k*k - k*l + l*l))
                        
                    elif i != 0 and k == 0:
                        for i_sign in range(2):
                            signed_i = i*(-1)**i_sign
                            m_eisenstein = eisenstein.EisensteinInt(j, signed_i)
                            n_eisenstein = eisenstein.EisensteinInt(l, k)
                            m_n_gcd = m_eisenstein.gcd(n_eisenstein)
                            _, m_arg = m_eisenstein.polar_form()
                            _, n_arg = n_eisenstein.polar_form()
                            arg_diff = m_arg - n_arg
                            if m_n_gcd.is_unit() and (i*i - signed_i*j + j*j > k*k - k*l + l*l or (i*i - signed_i*j + j*j == k*k - k*l + l*l and arg_diff > 0)):
                                ratio_dict[(j, signed_i, l, k)] = math.sqrt((i*i - signed_i*j + j*j)/(k*k - k*l + l*l))

                    elif i == 0 and k != 0:
                        for k_sign in range(2):
                            signed_k = k*(-1)**k_sign
                            m_eisenstein = eisenstein.EisensteinInt(j, i)
                            n_eisenstein = eisenstein.EisensteinInt(l, signed_k)
                            m_n_gcd = m_eisenstein.gcd(n_eisenstein)
                            _, m_arg = m_eisenstein.polar_form()
                            _, n_arg = n_eisenstein.polar_form()
                            arg_diff = m_arg - n_arg
                            if m_n_gcd.is_unit() and (i*i - i*j + j*j > k*k - signed_k*l + l*l or (i*i - i*j + j*j == k*k - signed_k*l + l*l and arg_diff > 0)):
                                ratio_dict[(j, i, l, signed_k)] = math.sqrt((i*i - i*j + j*j)/(k*k - signed_k*l + l*l))

                    elif i != 0 and k != 0:
                        for i_sign in range(2):
                            for k_sign in range(2):
                                signed_i = i*(-1)**i_sign
                                signed_k = k*(-1)**k_sign
                                m_eisenstein = eisenstein.EisensteinInt(j, signed_i)
                                n_eisenstein = eisenstein.EisensteinInt(l, signed_k)
                                _, m_arg = m_eisenstein.polar_form()
                                _, n_arg = n_eisenstein.polar_form()
                                arg_diff = m_arg - n_arg
                                m_n_gcd = m_eisenstein.gcd(n_eisenstein)
                                if m_n_gcd.is_unit() and (i*i - signed_i*j + j*j > k*k - signed_k*l + l*l or (i*i - signed_i*j + j*j == k*k - signed_k*l + l*l and arg_diff >= 0)): # = for (1, 0, 1, 0)
                                    ratio_dict[(j, signed_i, l, signed_k)] = math.sqrt((i*i - signed_i*j + j*j)/(k*k - signed_k*l + l*l))
                                
    sorted_ratio_dict = sorted(ratio_dict.items(), key=operator.itemgetter(1))
    return sorted_ratio_dict

def main():
    if len(sys.argv) == 3:
        filename = sys.argv[1]
        max_index = int(sys.argv[2])
    else:
        filename = 'eisenstein_list.py'
        max_index = 61
    if not os.path.isfile(filename):
        dict = make_eisenstein_supercell_dict(max_index)

        with open(filename, 'w') as f:
            for i, item in enumerate(dict):
                if i == 0:
                    f.write(f'sorted_dict = [{str(item)},\n')
                else:
                    f.write(f'{str(item)},\n')
            f.write(']\n')

def find_ratio(target_value, tol=1e-2):
    def make_eisenstein_matrix(m11, m12):
        return np.array([[m11, m12], [-m12, m11 - m12]])

    def make_return_dict(i):
        m = sorted_dict[i][0][0] + sorted_dict[i][0][1]*omega
        n = sorted_dict[i][0][2] + sorted_dict[i][0][3]*omega
        angle = cmath.phase(m/n)/math.pi*180

        return_dict = {'Transform_a': make_eisenstein_matrix(sorted_dict[i][0][0], sorted_dict[i][0][1]),
            'Transform_b': make_eisenstein_matrix(sorted_dict[i][0][2], sorted_dict[i][0][3]),
            'Quad_int_a': [sorted_dict[i][0][0], sorted_dict[i][0][1]],
            'Quad_int_b': [sorted_dict[i][0][2], sorted_dict[i][0][3]],
            'angle': angle,
            'table_ratio': sorted_dict[i][1]}
        return return_dict

    omega = -0.5+cmath.sqrt(-3)/2

    sorted_ratio_value = [d[1] for d in sorted_dict]
    n_left = bisect.bisect_left(sorted_ratio_value, target_value)

    i = n_left
    return_dict_list = []
    return_dict_list.append(make_return_dict(i))
    i = i + 1
    while (abs(sorted_ratio_value[n_left] - sorted_ratio_value[i]) < tol*target_value):
        return_dict_list.append(make_return_dict(i))
        i = i + 1
    
    return return_dict_list

def output_ratio_dict_list(ratio_dict_list, target_ratio):
    def str_ratio_dict(ratio_dict):
        abs_eisenstein_a = ratio_dict["Quad_int_a"][0]**2 - ratio_dict["Quad_int_a"][0]*ratio_dict["Quad_int_a"][1] + ratio_dict["Quad_int_a"][1]**2
        return f'({ratio_dict["Quad_int_a"][0]}+{ratio_dict["Quad_int_a"][1]}w)/({ratio_dict["Quad_int_b"][0]}+{ratio_dict["Quad_int_b"][1]}w), {ratio_dict["angle"]}, {abs_eisenstein_a}, {1 - ratio_dict["table_ratio"]/target_ratio}'

    output_str_list = []
    for dict_element in ratio_dict_list:
        output_str_list.append(str_ratio_dict(dict_element))

    return output_str_list

def output_file_ratio_dict_list(ratio_dict_list, target_ratio, filename):
    with open(filename, 'w') as f:
        f.writelines('\n'.join(output_ratio_dict_list(ratio_dict_list, target_ratio)))

def print_common_supercell_transform_dict(transform_dict, test_r):
    transform_dict_str = f'[{transform_dict["Transform_a"][0][0]: 03d} {transform_dict["Transform_a"][0][1]: 03d}]l_a = exp(i*{transform_dict["angle"]:0.02f}/180*pi)[{transform_dict["Transform_b"][0][0]: 03d} {transform_dict["Transform_b"][0][1]: 03d}]({transform_dict["table_ratio"]/test_r})l_b ({transform_dict["table_ratio"]:.4f})\n[{transform_dict["Transform_a"][1][0]: 03d} {transform_dict["Transform_a"][1][1]: 03d}]                            [{transform_dict["Transform_b"][0][0]: 03d} {transform_dict["Transform_b"][1][1]: 03d}]'
    print(transform_dict_str) 
    m1 = np.array(transform_dict["Transform_a"])
    m2 = np.array(transform_dict["Transform_b"])
    
    #omega = -0.5+cmath.sqrt(-3)/2
    #v1 = np.array([[1.0], [omega]])
    #print(m1@v1)
    #print(test_r*cmath.exp(1j*transform_dict["angle"]/180.0*math.pi)*m2@v1)
    #print(transform_dict['table_ratio']*cmath.exp(1j*transform_dict["angle"]/180.0*math.pi)*m2@v1)

if __name__ == '__main__':
    main()
