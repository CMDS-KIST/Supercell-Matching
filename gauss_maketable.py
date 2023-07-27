import math, cmath, bisect
import operator
import quadraticinteger

import numpy as np
import os
import sys

def make_gaussian_supercell_dict(max_index):
    #omega = cmath.sqrt(-1)
    ratio_dict = {}
    for i in range(0, max_index):
        for j in range(i+1, max_index):
            for k in range(0, max_index):
                for l in range(0, max_index): # if k > l, we miss (2, -1, 1, -2). 2 - i and 1 - 2i are coprime
                    if i == 0 and k == 0:
                        if l != 0:
                            m_gaussian = quadraticinteger.GaussInt_square(j, i)
                            n_gaussian = quadraticinteger.GaussInt_square(l, k)
                            m_n_gcd = m_gaussian.gcd(n_gaussian)
                            m_arg = m_gaussian.arg()
                            n_arg = n_gaussian.arg()
                            arg_diff = m_arg - n_arg
                            if m_n_gcd.is_unit() and (i*i + j*j > k*k + l*l or (i*i + j*j == k*k + l*l and arg_diff >= 0)):
                                ratio_dict[(j, i, l, k)] = math.sqrt((i*i + j*j)/(k*k + l*l))
                        
                    elif i != 0 and k == 0:
                        for i_sign in range(2):
                            signed_i = i*(-1)**i_sign
                            m_gaussian = quadraticinteger.GaussInt_square(j, signed_i)
                            n_gaussian = quadraticinteger.GaussInt_square(l, k)
                            m_n_gcd = m_gaussian.gcd(n_gaussian)
                            m_arg = m_gaussian.arg()
                            n_arg = n_gaussian.arg()
                            arg_diff = m_arg - n_arg
                            if m_n_gcd.is_unit() and (i*i + j*j > k*k + l*l or (i*i + j*j == k*k + l*l and arg_diff > 0)):
                                ratio_dict[(j, signed_i, l, k)] = math.sqrt((i*i + j*j)/(k*k + l*l))

                    elif i == 0 and k != 0:
                        for k_sign in range(2):
                            signed_k = k*(-1)**k_sign
                            m_gaussian = quadraticinteger.GaussInt_square(j, i)
                            n_gaussian = quadraticinteger.GaussInt_square(l, signed_k)
                            m_n_gcd = m_gaussian.gcd(n_gaussian)
                            m_arg = m_gaussian.arg()
                            n_arg = n_gaussian.arg()
                            arg_diff = m_arg - n_arg
                            if m_n_gcd.is_unit() and (i*i + j*j > k*k + l*l or (i*i + j*j == k*k + l*l and arg_diff > 0)):
                                ratio_dict[(j, i, l, signed_k)] = math.sqrt((i*i + j*j)/(k*k + l*l))

                    elif i != 0 and k != 0:
                        for i_sign in range(2):
                            for k_sign in range(2):
                                signed_i = i*(-1)**i_sign
                                signed_k = k*(-1)**k_sign
                                m_gaussian = quadraticinteger.GaussInt_square(j, signed_i)
                                n_gaussian = quadraticinteger.GaussInt_square(l, signed_k)
                                m_arg = m_gaussian.arg()
                                n_arg = n_gaussian.arg()
                                arg_diff = m_arg - n_arg
                                m_n_gcd = m_gaussian.gcd(n_gaussian)
                                if m_n_gcd.is_unit() and (i*i + j*j > k*k + l*l or (i*i + j*j == k*k + l*l and arg_diff >= 0)): # = for (1, 0, 1, 0)
                                    ratio_dict[(j, signed_i, l, signed_k)] = math.sqrt((i*i + j*j)/(k*k + l*l))
                                
    sorted_ratio_dict = sorted(ratio_dict.items(), key=operator.itemgetter(1))
    return sorted_ratio_dict

def main():
    if len(sys.argv) == 3:
        filename = sys.argv[1]
        max_index = int(sys.argv[2])
    else:
        filename = 'gaussian_list_pos_arg.py'
        max_index = 61
    if not os.path.isfile(filename):
        dict = make_gaussian_supercell_dict(max_index)

        with open(filename, 'w') as f:
            for i, item in enumerate(dict):
                if i == 0:
                    f.write(f'sorted_dict = [{str(item)},\n')
                else:
                    f.write(f'{str(item)},\n')
            f.write(']\n')

if __name__ == "__main__":
    main()