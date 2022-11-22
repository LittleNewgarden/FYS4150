import numpy as np

# lists of times taken from terminal
time_parallel = np.array([0.195661, 0.192198, 0.191404, 0.20337, 0.198426, 0.192719, 0.198636, 0.192521, 0.20017, 0.195302])
time_non_parallel = np.array([1.74136, 1.77473, 1.73817, 1.71437, 1.78727, 1.75898, 1.72684, 1.77351, 1.74577, 1.76404])

#averages
par =  np.sum(time_parallel)/len(time_parallel)
non_par = np.sum(time_non_parallel)/len(time_non_parallel)

print("average time parrallel : ", par)
print("average time non-parrallel : ", non_par)

print("ratio : ", par/non_par)