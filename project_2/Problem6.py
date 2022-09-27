import numpy as np
import matplotlib.pyplot as plt


plt.style.use(["bmh"])
plt.rcParams.update({"font.size": 11})#,


# oppg 6 a)

x = np.linspace(0,1,11)

v = np.zeros((3,11))
a = np.zeros((3,11))

v = np.loadtxt("prob6_11.txt",unpack=True)

a = np.loadtxt("prob6_analytical_11.txt",unpack=True)



for i in range(len(v)):
    plt.plot(x,v[i])
    plt.plot(x,a[i])

    plt.legend(["discretized n_steps=10","analytical n_steps=10"])
    plt.ylabel(r'$v(\hat{x})/u(\hat{x})$')
    plt.xlabel(r'$\hat{x}$')
    plt.savefig(f'figur_task_6_a_{i}.pdf')
    plt.show()










# oppg 6 b)


x = np.linspace(0,1,101)

v = np.zeros((3,101))
a = np.zeros((3,101))

v = np.loadtxt("prob6_101.txt",unpack=True)

a = np.loadtxt("prob6_analytical_101.txt",unpack=True)



for i in range(len(v)):
    plt.plot(x,v[i])
    plt.plot(x,a[i])

    plt.legend(["discretized n_steps=100","analytical n_steps=100"])
    plt.ylabel(r'$v(\hat{x})/u(\hat{x})$')
    plt.xlabel(r'$\hat{x}$')
    plt.savefig(f'figur_task_6_b_{i}.pdf')
    plt.show()



