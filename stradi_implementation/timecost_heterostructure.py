import time
import heterostructure

def get_function_timecost(func, args):
    begin_time = time.time()
    func(*args)
    end_time = time.time()
    return end_time - begin_time

time_index_list = [[], []]
max_index = 5
index_list = list(range(1, max_index))

for i in range(1, max_index):
    heterostructure_builder = heterostructure.heterostructure(1.24, 1.24, 90, 1, 1, 90, max_index, 1)

    time_trial_list = []
    func_list = [heterostructure_builder.find_quadratic_heterostructures]
    for j in range(1):
        method_time_trial_list = []
        for i in range(1):
            method_time_trial_list.append(get_function_timecost(func_list[j], [0.01]))
        time_trial_list.append(method_time_trial_list)

        time_index_list[j].append(sum(time_trial_list[j])/len(time_trial_list[j]))

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2)
for j in range(1):
    axes[j].plot(index_list, time_index_list[j])

plt.show()

