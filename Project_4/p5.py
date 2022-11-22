import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

quantity_support()
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 8)) # Larger figure sizes
plt.rc('font', size=18)

#     ---------- PROBLEM 5 ----------



names = ["20x20_per=50_n=10^6_all.txt", "20x20_per=100_n=10^6_all.txt"]

N = 20*20  # number of spins

start = 0
end = 3*10**4  #area of intrest in the files

# extracting data
EM_per_50 = np.loadtxt(names[0])   
EM_per_100 = np.loadtxt(names[1])

E_50 = EM_per_50[start:end,0::2]/N

M_50 = EM_per_50[start:end,1::2]/N

E_100 = EM_per_100[start:end,0::2]/N
M_100 = EM_per_100[start:end,1::2]/N

NN = end-start


avrg_e_50, avrg_m_50 = np.zeros((NN,2)), np.zeros((NN,2))
avrg_e_100, avrg_m_100 = np.zeros((NN,2)), np.zeros((NN,2))

avrg_e_50[0], avrg_m_50[0] = E_50[0], np.abs(M_50[0])
avrg_e_100[0], avrg_m_100[0] =  E_100[0], np.abs(M_100[0])


# findign the cumalative averages
for i in range(NN-1):
    avrg_e_50[i+1] = ((i+1)*avrg_e_50[i] + E_50[i+1])/(i+2)
    avrg_m_50 [i+1]= ((i+1)*avrg_m_50[i] + np.abs(M_50[i+1]))/(i+2)

    avrg_e_100[i+1] =  ((i+1)*avrg_e_100[i] + E_100[i+1])/(i+2)
    avrg_m_100[i+1] =  ((i+1)*avrg_m_100[i] + np.abs(M_100[i+1]))/(i+2)

NN_range= np.arange(NN)




plt.plot(NN_range, avrg_e_100)#, label="100%")
plt.plot(NN_range, avrg_e_50)#, label="50%")
plt.legend(["T=1, 100%", "T=2.4, 100%", "T=1, 50%", "T=2.4, 50%"])
plt.xlabel("Monte Carlo cycles")
plt.ylabel(r"$\langle \varepsilon \rangle$ [J/N]")
plt.savefig("figs/p5_e.pdf")
plt.show()



plt.plot(NN_range, avrg_m_100)#, label="100%")
plt.plot(NN_range, avrg_m_50)#, label="50%")
plt.legend(["T=1, 100%", "T=2.4, 100%", "T=1, 50%", "T=2.4, 50%"])
plt.xlabel("Monte Carlo cycles")
plt.ylabel(r"$\langle |m| \rangle$ [1/N]")
plt.savefig("figs/p5_m.pdf")
plt.show()







