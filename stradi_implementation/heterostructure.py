import numpy as np
import math
import cmath
import bisect
import time

from sympy import fu

class KeyList(object):
    def __init__(self, l, key):
        self.l = l
        self.key = key

    def __len__(self):
        return len(self.l)

    def __getitem__(self, index):
        return self.key(self.l[index])

def rotation_matrix(angle):
    return np.array([[math.cos(angle/180.0*math.pi), math.sin(-angle/180.0*math.pi)], [math.sin(angle/180.0*math.pi), math.cos(angle/180.0*math.pi)]])

def quadratic_integer_matrix(m, n, omega):
    return np.array([[m, 0], [0, m]]) + np.array([[0, n], [-n*omega*omega.conjugate(), n*(omega + omega.conjugate())]])

class heterostructure:
    def __init__(self, top_a, top_b, top_gamma, bottom_a, bottom_b, bottom_gamma, max_index, angle_step):
        self.top_a = top_a
        self.top_b = top_b
        self.top_gamma = top_gamma
        self.bottom_a = bottom_a
        self.bottom_b = bottom_b
        self.bottom_gamma = bottom_gamma
        self.max_index = max_index
        self.angle_step = angle_step
        self.a_top = np.array([[self.top_a], [0.0]])
        self.b_top = np.array([[self.top_b*np.cos(top_gamma/180.0*math.pi)], [self.top_b*np.sin(top_gamma/180.0*math.pi)]])
        self.lattice_matrix_top = (np.concatenate([self.a_top, self.b_top], 1))
        self.a_bottom = np.array([[self.bottom_a], [0.0]])
        self.b_bottom = np.array([[self.bottom_b*np.cos(bottom_gamma/180.0*math.pi)], [self.bottom_b*np.sin(bottom_gamma/180.0*math.pi)]])
        self.e1_bottom = self.a_bottom/(np.linalg.norm(self.a_bottom)**2)
        self.e2_bottom = self.b_bottom/(np.linalg.norm(self.b_bottom)**2)
        self.lattice_matrix_bottom = (np.concatenate([self.a_bottom, self.b_bottom], 1))
        self.e_matrix_bottom = (np.concatenate([self.e1_bottom, self.e2_bottom], 1))
        self.quadratic_ratio_table = None

    def find_heterostructures(self, strain_max = 1.0):
        supercell_list = []
        for angle in range(0, 360, self.angle_step):
            rotation_matrix = np.array([[math.cos(angle/180.0*math.pi), math.sin(-angle/180.0*math.pi)], [math.sin(angle/180.0*math.pi), math.cos(angle/180.0*math.pi)]])
            for i in range(-self.max_index, self.max_index + 1):
                for j in range(0, self.max_index + 1):
                    for k in range(-self.max_index, self.max_index + 1):
                        for l in range(1, self.max_index + 1):
                            T1 = np.array([[i, j], [k, l]])
                            if abs(np.linalg.det(T1)) > 1e-4:
                                lattice_matrix_supercell_top = rotation_matrix@(self.lattice_matrix_top@np.transpose(T1))
                                coefficient_matix_lattice_matrix_supercell_top_projected_to_bottom = np.transpose(lattice_matrix_supercell_top)@self.e_matrix_bottom
                                T2 = np.round(coefficient_matix_lattice_matrix_supercell_top_projected_to_bottom)
                                lattice_matrix_supercell_bottom = self.lattice_matrix_bottom@np.transpose(T2)
                                e11 = abs(lattice_matrix_supercell_top[0, 0]/lattice_matrix_supercell_bottom[0, 0]) - 1
                                e22 = abs(lattice_matrix_supercell_top[1, 1]/lattice_matrix_supercell_bottom[1, 1]) - 1
                                e12 = (lattice_matrix_supercell_top[1, 0] - lattice_matrix_supercell_top[0, 0]/lattice_matrix_supercell_bottom[0, 0]*lattice_matrix_supercell_bottom[1, 0])/2.0/lattice_matrix_supercell_bottom[1, 1]
                                #strain_tensor = np.array([e11, e12], [e12, e22])
                                abs_stress = (abs(e11) + abs(e22) + abs(e12))/3.0
                                if abs_stress < strain_max:
                                    supercell_list.append([T1, T2, angle, abs_stress])
        return supercell_list

    def construct_quadratic_table(self):
        omega_top = self.top_a/self.top_b*cmath.exp(1j*self.top_gamma)
        omega_bottom = self.bottom_b/self.bottom_b*cmath.exp(1j*self.bottom_gamma)
        if abs(omega_bottom - omega_top) > 1e-4:
            exit(1)
        else:
            omega = omega_top

        self.quadratic_ratio_table = []

        for i in range(-self.max_index, self.max_index + 1):
            for j in range(0, self.max_index + 1):
                for k in range(-self.max_index, self.max_index + 1):
                    for l in range(1, self.max_index + 1):
                        if not ((i == 0 and j == 0) or (k == 0 and l == 0)):
                            gen_eigenvalue_T1T2 = (k + l*omega)/(i + j*omega)
                            r = abs(gen_eigenvalue_T1T2)
                            angle = cmath.phase(gen_eigenvalue_T1T2)
                            self.quadratic_ratio_table.append([i, j, k, l, r, angle])

        sorted(self.quadratic_ratio_table, key = lambda x: x[5])
        sorted(self.quadratic_ratio_table, key = lambda x: x[4])
    
    def find_quadratic_heterostructures(self, strain_max = 1.0):
        omega_top = self.top_a/self.top_b*cmath.exp(1j*self.top_gamma)
        omega_bottom = self.bottom_b/self.bottom_b*cmath.exp(1j*self.bottom_gamma)
        if abs(omega_bottom - omega_top) > 1e-4:
            exit(1)
        else:
            omega = omega_top

        supercell_list = []
        if not self.quadratic_ratio_table:
            quadratic_ratio_table = []
            #s_time = time.time()
            for i in range(-self.max_index, self.max_index + 1):
                for j in range(0, self.max_index + 1):
                    for k in range(-self.max_index, self.max_index + 1):
                        for l in range(1, self.max_index + 1):
                            if not ((i == 0 and j == 0) or (k == 0 and l == 0)):
                                gen_eigenvalue_T1T2 = (k + l*omega)/(i + j*omega)
                                r = abs(gen_eigenvalue_T1T2)
                                angle = cmath.phase(gen_eigenvalue_T1T2)
                                quadratic_ratio_table.append([i, j, k, l, r, angle])

            sorted(quadratic_ratio_table, key = lambda x: x[5])
            sorted(quadratic_ratio_table, key = lambda x: x[4])
        else:
            quadratic_ratio_table = self.quadratic_ratio_table

        target_ratio = self.bottom_a/self.top_a
        #e_time = time.time()
        #print(f'Table time: {e_time - s_time:.3f}')
        
        import operator
        #closest_ratio_item_index = bisect.bisect_left(quadratic_ratio_table, target_ratio, key = lambda x: x[4])
        #s_time = time.time()
        for i in range(1000):
            closest_ratio_item_index = bisect.bisect_left(KeyList(quadratic_ratio_table, operator.itemgetter(4)), target_ratio)
        closest_ratio_item = quadratic_ratio_table[closest_ratio_item_index]
        stress = target_ratio - closest_ratio_item[4]
        #e_time = time.time()
        #print(f'Find time: {e_time - s_time:.3e}')

        #s_time = time.time()
        if stress < strain_max:
            supercell_list.append([quadratic_integer_matrix(closest_ratio_item[0], closest_ratio_item[1], omega), quadratic_integer_matrix(closest_ratio_item[2], closest_ratio_item[3], omega), closest_ratio_item[5], stress])

        ratio_index = closest_ratio_item_index + 1
        while (True):
            next_item = quadratic_ratio_table[ratio_index]
            stress = target_ratio - next_item[4]
            if stress < strain_max:
                supercell_list.append([quadratic_integer_matrix(next_item[0], next_item[1], omega), quadratic_integer_matrix(next_item[2], next_item[3], omega), next_item[5], stress])
                ratio_index = ratio_index + 1
            else:
                break
        #e_time = time.time()
        #print(f'Find time: {e_time - s_time:.3e}, number of solutions = {len(supercell_list)}')
        return supercell_list

def get_function_timecost(func, args):
    begin_time = time.time()
    func(*args)
    end_time = time.time()
    return end_time - begin_time

if __name__ == '__main__':
    time_q_list = []
    time_f_list = []
    max_index = 10
    range_for_index_scan = range(1, max_index + 1)
    list_for_index_scan = list(range_for_index_scan)
    for i in range_for_index_scan:
        test_object = heterostructure(1.24, 1.24, 90, 1, 1, 90, i, 1)
        test_object.construct_quadratic_table()
        
        time_q_list.append(get_function_timecost(test_object.find_quadratic_heterostructures, [0.01])/1000)
    #for i in range(1, 11):
    #    test_object = heterostructure(1.24, 1.24, 90, 1, 1, 90, i, 1)
    #    time_f_list.append(get_function_timecost(test_object.find_heterostructures, [0.01]))

import matplotlib.pyplot as plt
from functools import partial

fig, axes = plt.subplots(1, 2)
axes[0].plot(list_for_index_scan, time_q_list, label='Quadratic integer table')
axes[1].plot(list_for_index_scan, time_q_list, label='Quadratic integer table')
#plt.plot(list(range(1, 11)), time_f_list, label='Scanning')

axes[1].set_yscale('function', functions=(partial(np.power, 10.0), np.log10))
plt.show()

with open('time_table_notablegen.txt', 'w') as f:
    f.write('Index, Quadratic_integer_table\n')
    for i in range(1, 11):
        f.write(f'{i}, {time_q_list[i-1]}\n')
