import numpy as np
import matplotlib.pyplot as plt

#oepn and read text file containing x and u(x) values
with open("u_x_output.txt", "r+") as file:
    words = []
    for line in file:
        words.append([float(line.split()[0]), float(line.split()[1])])
        
    words = (np.array(words))



#plotting u(x) as a function of x

plt.plot(words[:,0],words[:,1])
plt.xlabel("x")
plt.ylabel("u(x)")
plt.savefig("figure_task_2.pdf")
plt.show()