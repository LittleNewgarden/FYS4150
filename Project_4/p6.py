import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

quantity_support()
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 8)) # Larger figure sizes
plt.rc('font', size=18)

#     ---------- PROBLEM 6 ----------



names = ["20x20_per=50_n=10^6_all.txt"]


N = 20*20

start = 2*10**4  # removing burn-in
end = -1


EM_per_50 = np.loadtxt(names[0])

e = EM_per_50[start:end,0::2]/N

maxdE = 4
maxde = maxdE/N #bin width

nbins = np.array([int(np.abs(np.max(e[:,0]) - np.min(e[:,0]) )/maxde), int(np.abs(np.max(e[:,1]) - np.min(e[:,1]) )/maxde)])

bin_edges0 = np.linspace(np.min(e[:,0]) - maxde/2, np.max(e[:,0]) + maxde/2,nbins[0]+2)
bin_edges1 = np.linspace(np.min(e[:,1]) - maxde/2, np.max(e[:,1]) + maxde/2,nbins[1]+2)


plt.hist(e[:,0],bins=bin_edges0,density=True)
plt.xlabel(r"$\varepsilon $ [J/N]")
plt.ylabel("Column height")
plt.savefig("figs/p6_T=1.pdf")
plt.show()


plt.hist(e[:,1],bins=bin_edges1,density=True)
plt.xlabel(r"$\varepsilon $ [J/N]")
plt.ylabel("Column height")
plt.savefig("figs/p6_T=2.4.pdf")
plt.show() 