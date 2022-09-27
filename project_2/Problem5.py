import numpy as np
import matplotlib.pyplot as plt


plt.style.use(["bmh"])
plt.rcParams.update({"font.size": 15})#,


# reading file
N, iterations = np.loadtxt("iterations_per_N.txt",unpack=True)


# plotting size of matirx versus the number of iterations needed
plt.plot(N,iterations)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('iterations')
plt.savefig(f'figur_task_5.pdf')
plt.show()