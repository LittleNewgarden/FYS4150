import numpy as np
import matplotlib.pyplot as plt


with open("u_x_output.txt", "r+") as file:
    words = []
    for line in file:
        words.append([float(line.split()[0]), float(line.split()[1])])
        
    
    words = (np.array(words))
    print(words)

plt.plot(words[:,0],words[:,1])

plt.show()