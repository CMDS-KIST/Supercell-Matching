import math, cmath
from quadraticinteger import EuclideanDomainInt, DedekindDomainInt
import pickle
import argparse

quadratic_field_norm_euclidean = [-11, -7, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73]

def make_quadratic_supercell_dict(max_index, d):
    omega = cmath.sqrt(-1)
    
    ratio_dict = {}

    d_mod4 = divmod(d, 4)[1]
    if d_mod4 == 3 or d_mod4 == 2:
        omega = cmath.sqrt(d)
    elif d_mod4 == 1:
        omega = (-1 + cmath.sqrt(d))/2
    
    if d == -1 or d == -3:
        min_index = 0 # Eisenstein integers and Gaussian integers have more than 2 units
    else:
        min_index = -max_index # The other quadratic integers have only 2 units, 1 and -1

    normsquare_func = lambda x, y: ((x + omega*y)*(x + omega.conjugate()*y)).real

    if d in quadratic_field_norm_euclidean and d != -3:
        for i in range(1, max_index):
            for j in range(min_index, max_index):
                for k in range(1, max_index):
                    for l in range(min_index, max_index):
                        normsquare_ij = normsquare_func(i, j)
                        normsquare_kl = normsquare_func(k, l)
                        z_m = EuclideanDomainInt(i, j, d)
                        z_n = EuclideanDomainInt(k, l, d)
                        rotation_angle = (cmath.phase(z_m.complex_form()) - cmath.phase(z_n.complex_form()))/math.pi*180
                        if z_m.gcd(z_n).is_unit() and not ((i == 0 and j == 0) or (k == 0 and l == 0)) and normsquare_ij >= normsquare_kl:
                            ratio_dict[(i, j, k, l)] = [math.sqrt(normsquare_ij/normsquare_kl), rotation_angle]
    elif d == -3:
        for i in range(1, max_index):
            for j in range(min_index, max_index):
                for k in range(1, max_index):
                    for l in range(min_index, max_index):
                        normsquare_ij = normsquare_func(i + j, j) # Primes of eisenstein integers on first sextant (i + j*(1 + \omega))
                        normsquare_kl = normsquare_func(k + l, l)
                        z_m = EuclideanDomainInt(i + j, j, d)
                        z_n = EuclideanDomainInt(k + l, l, d)
                        rotation_angle = (cmath.phase(z_m.complex_form()) - cmath.phase(z_n.complex_form()))/math.pi*180
                        if z_m.gcd(z_n).is_unit() and not ((i == 0 and j == 0) or (k == 0 and l == 0)) and normsquare_ij >= normsquare_kl:
                            ratio_dict[(i + j, j, k + l, l)] = [math.sqrt(normsquare_ij/normsquare_kl), rotation_angle]
    else:
        for i in range(1, max_index):
            for j in range(min_index, max_index):
                for k in range(1, max_index):
                    for l in range(min_index, max_index):
                        normsquare_ij = normsquare_func(i, j) # Primes of eisenstein integers on first sextant (i + j*(1 + \omega))
                        normsquare_kl = normsquare_func(k, l)
                        z_m = DedekindDomainInt(i, j, d)
                        z_n = DedekindDomainInt(k, l, d)
                        rotation_angle = (cmath.phase(z_m.complex_form()) - cmath.phase(z_n.complex_form()))/math.pi*180
                        if z_m.gcd(z_n).is_unit() and not ((i == 0 and j == 0) or (k == 0 and l == 0)) and normsquare_ij >= normsquare_kl:
                            ratio_dict[(i, j, k, l)] = [math.sqrt(normsquare_ij/normsquare_kl), rotation_angle]

    sorted_ratio_dict_1 = sorted(ratio_dict.items(), key=lambda item: normsquare_func(item[0][0], item[0][1]))
    sorted_ratio_dict = sorted(sorted_ratio_dict_1, key=lambda item: item[1][0])
    return sorted_ratio_dict

def dump_table(sorted_ratio_dict, filename):
    with open(filename, 'wb') as f:
        pickle.dump(sorted_ratio_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

def write_table_as_list(sorted_ratio_list, filename):
    with open(filename, 'w') as f:
            for j, item in enumerate(sorted_ratio_list):
                if j == 0:
                    f.write(f'sorted_ratio_dict = [{str(item)},\n')
                else:
                    f.write(f'{str(item)},\n')
            
            f.write(']\n')

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(prog='Zquad_maketable.py', description='Make a table to find common supercells for a quadatic integer lattice')
    argparser.add_argument('D', type=int, help='D for a quadratic integer O[Z[D]] (for two-dimensional lattice, D should be less than 0)')
    argparser.add_argument('--max_index', type=int, default=60, help='Maximum limit of the coefficients of the quadratic integers')

    args = argparser.parse_args()
    max_index = args.max_index
    quad_d = args.D
    filename_pickle = f'OZsqrt{-quad_d}i_list_{max_index}.pickle'

    if quad_d < 0:
        sorted_ratio_list = make_quadratic_supercell_dict(max_index, quad_d)
    
        dump_table(sorted_ratio_list, filename_pickle)
        print(f'A table file {filename_pickle} was written.')
    else:
        print('D should be less than 0.')
        exit(1)