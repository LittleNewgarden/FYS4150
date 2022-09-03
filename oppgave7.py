import numpy as np
import matplotlib.pyplot as plt


# with open("u_x_output.txt", "r+") as file:
#     words = []
#     for line in file:
#         words.append([float(line.split()[0]), float(line.split()[1])])
        
    
#     words = (np.array(words))
#     print(words)

# plt.plot(words[:,0],words[:,1])

# plt.show()


def u(x):
    return 1- (1 - np.exp(-10))*x - np.exp(-10*x)

x = []
v = []
n = [10,100,1000,100000]
for i in n:
    x.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[0])
    v.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[1])

plt.plot(x[0],v[0],label='10 steps')
plt.plot(x[1],v[1],label='100 steps')
plt.plot(x[2],v[2],label='1000 steps')
plt.plot(x[3],v[3],label='100000 steps')
plt.plot(x[3],u(x[3]), label='u(x)')
plt.legend()
plt.show()