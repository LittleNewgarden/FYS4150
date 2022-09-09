import numpy as np
import matplotlib.pyplot as plt


# reading the text file produced by oppgave10.cpp
general, special =  (np.loadtxt("general_vs_special.txt",unpack=True))


n_steps = [10,10**2,10**3,10**4,10**5,10**6]

# plotting the time used by the general and special algorithem for different n_steps

plt.plot(n_steps,general, label="general")
plt.plot(n_steps,special, label="special")
plt.xlabel("n_steps")
plt.ylabel("time used")
plt.legend()
plt.xscale("log")
plt.savefig('figur_task_10.pdf')
plt.show()