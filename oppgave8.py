import numpy as np
import matplotlib.pyplot as plt

# A)

# exact solution of u(x)
def u(x):
    return 1- (1 - np.exp(-10))*x - np.exp(-10*x)

# reading all the files produced by oppgave7.cpp and extracting x and v
x = []
v = []
n = [10**1,10**2,10**3,10**4,10**5,10**6,10**7]
for i in n:
    x.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[0])
    v.append(np.loadtxt(f"P7_{i}.txt",unpack=True)[1])


# calculating the absolute error using n=10 to n = 10^4
abs_error = []
N = len(n)
for i in range(N-3):
    U = u(x[i][1:-1]) #calculating the exact solution at the points x_i
    abs_error.append( np.abs(U - v[i][1:-1])) 


# plotting th absolute error for n up to 10^4

plt.title("absolute error" )
plt.plot(x[0][1:-1],abs_error[0],label='10 steps')
plt.plot(x[1][1:-1],abs_error[1],label='100 steps')
plt.plot(x[2][1:-1],abs_error[2],label='1000 steps')
plt.plot(x[3][1:-1],abs_error[3],label='10000 steps')
plt.yscale("log")
plt.xlabel('x')
plt.ylabel('absolute error')
plt.legend()
plt.savefig('figur_task_8_a.pdf')
plt.show()

# B)

# calculating the absolute error using n up to 10^7
rel_error = []
for i in range(N):
    U = u(x[i][1:-1]) #calculating the exact solution at the points x_i
    rel_error.append( np.abs((U - v[i][1:-1])/U))


# plotting th absolute error for n up to 10^4

plt.title("relative error" )
plt.plot(x[0][1:-1],rel_error[0],label='10 steps')
plt.plot(x[1][1:-1],rel_error[1],label='100 steps')
plt.plot(x[2][1:-1],rel_error[2],label='1000 steps')
plt.plot(x[3][1:-1],rel_error[3],label='10000 steps')
plt.yscale("log")
plt.xlabel('x')
plt.ylabel('relative error')
plt.legend()
plt.savefig('figur_task_8_b.pdf')
plt.show()


# C)


# making a table of the max value of the relative errors for each n
max_rel = []
print("n_steps    relative error")
for i in range(len(n)):
    max_rel.append(np.max(rel_error[i]))
    print(n[i],"       ",np.max(rel_error[i]))


# plotting the max relative error as a function of n

plt.plot(n,max_rel,"o")
plt.xlabel("n_steps")
plt.ylabel("max relative error")
plt.yscale('log')
plt.xscale('log')
plt.savefig('figur_task_8_c.pdf')
plt.show()

