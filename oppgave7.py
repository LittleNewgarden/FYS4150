import numpy as np
import matplotlib.pyplot as plt


# exact solution of u(x)
def u(x):
    return 1- (1 - np.exp(-10))*x - np.exp(-10*x)

# reading files produced in oppgave7.cpp, only taking up to n=10^4
x = []
v = []
n = [10,100,1000,10000]
for i in n:
    x.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[0])
    v.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[1])

# plotting the approximation to u(x) together with u(x)

plt.plot(x[0],v[0],label='10 steps')
plt.plot(x[1],v[1],label='100 steps')
plt.plot(x[2],v[2],label='1000 steps')
plt.plot(x[3],v[3],label='10000 steps')
plt.plot(x[3],u(x[3]), label='u(x)')
plt.legend()
plt.savefig("figure_task_7.pdf")
plt.xlabel('x')
plt.ylabel('v(x) and g(x)')
plt.show()