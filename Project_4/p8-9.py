import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

quantity_support()
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 8)) # Larger figure sizes
plt.rc('font', size=18)



#     ---------- PROBLEM 8 ----------



names = ["40x40_n=10^6_Tn=100_average.txt", "60x60_n=10^6_Tn=100_average.txt", "80x80_n=10^6_Tn=100_average.txt", "100x100_n=10^6_Tn=100_average.txt"]

L = [40, 60, 80, 100] # different sizes

Tn = 100 # number of temperature steps

T_range = np.linspace(2.1,2.4,Tn) #temperature range

values = []
for j in range(len(names)):                # loading all values
    values.append(np.loadtxt(names[j]))

for i in range(len(names)):    # epsilon plots for different L
    plt.plot(T_range, values[i][:,0], label=f"L={L[i]}")
plt.xlabel(r"Temperature $[J/k_B]$")
plt.ylabel(r"$\langle \varepsilon \rangle$ [J/N]")
plt.legend()
plt.savefig("figs/p8_e.pdf")
plt.show()

for i in range(len(names)):   # m plots for different L
    plt.plot(T_range, values[i][:,1], label=f"L={L[i]}")
plt.xlabel(r"Temperature $[J/k_B]$")
plt.ylabel(r"$\langle |m| \rangle$ [1/N]")
plt.legend()
plt.savefig("figs/p8_m.pdf")
plt.show()

for i in range(len(names)):   # C_vV plots for different L
    plt.plot(T_range, values[i][:,4], label=f"L={L[i]}")
plt.xlabel(r"Temperature $[J/k_B]$")
plt.ylabel(r"$C_V$ $[k_B/N]$")
plt.legend()
plt.savefig("figs/p8_Cv.pdf")
plt.show()

for i in range(len(names)):   # chi plots for different L
    plt.plot(T_range, values[i][:,5], label=f"L={L[i]}")
plt.xlabel(r"Temperature $[J/k_B]$")
plt.ylabel(r"$\chi$ [1/(JN)]")
plt.legend()
plt.savefig("figs/p8_chi.pdf")
plt.show()



# C_V and chi values for different L
Cv =  np.array([values[0][:,4] , values[1][:,4], values[2][:,4], values[3][:,4]])
Xi =  np.array([values[0][:,5] , values[1][:,5], values[2][:,5], values[3][:,5]])

# j = np.argmin(np.abs(Cv[:,:] - np.max(Cv[:,:],axis=1)[:,np.newaxis]),axis=1)


# print(" L :    |    40     |    60    |    80    |    100    |")
# print("Tc(Cv) : ", T_range[j])

# j1 = np.argmin(np.abs(Xi[:,:] - np.max(Xi[:,:],axis=1)[:,np.newaxis]),axis=1)


# print(" L :    |    40     |    60    |    80    |    100    |")
# print("Tc(Xi) : ", T_range[j1])





from scipy.optimize import curve_fit



#  gaussian func used in fitting
def f(x,a,b,c):
    return a*np.exp(-(x-b)**2/c**2)




# indicies of T used for fitting the peaks
A = 40
B = 85

A1 = 45
B1 = 75

A2 = 50
B2 = 70

A3 = 50
B3 = 70

T40 = T_range[A:B]
T60 = T_range[A1:B1]
T80 = T_range[A2:B2]
T100 = T_range[A3:B3]

# fitting the curves
a,b = curve_fit(f,T40, Cv[0,A:B])
a1,b1 = curve_fit(f,T60, Cv[1,A1:B1])
a2,b2 = curve_fit(f,T80, Cv[2,A2:B2])
a3,b3 = curve_fit(f,T100, Cv[3,A3:B3])


# The fitted curves
Cv_40 = f(T40,a[0],a[1],a[2])
Cv_60 = f(T60,a1[0],a1[1],a1[2])
Cv_80 = f(T80,a2[0],a2[1],a2[2])
Cv_100 = f(T100,a3[0],a3[1],a3[2])


# plotting the fitted curves with the data
plt.plot(T40, Cv_40, color="blue", label=f"L={L[0]}: Gaussian fit")#T_range[A:B], Cv_40)
plt.plot(T_range, values[0][:,4], "o", color="blue", label=f"L={L[0]}: Data")
plt.plot(T60, Cv_60, color="orange", label=f"L={L[1]}: Gaussian fit")
plt.plot(T_range, values[1][:,4], "o", color="orange", label=f"L={L[1]}: Data")
plt.plot(T80, Cv_80, color="red", label=f"L={L[2]}: Gaussian fit")
plt.plot(T_range, values[2][:,4], "o", color="red", label=f"L={L[2]}: Data")
plt.plot(T100, Cv_100, color="black", label=f"L={L[3]}: Gaussian fit")
plt.plot(T_range, values[3][:,4], "o", color="black", label=f"L={L[3]}: Data")
plt.xlabel(r"Temperature $[J/k_{B}]$ ")
plt.ylabel(r"$C_V$ $[k_B/N]$")
plt.legend()
plt.savefig("figs/p8_Cv_fit.pdf")
plt.show()



# finding Tc for differnt L
j40 = np.argmin(np.abs( Cv_40 - np.max(Cv_40)  ))
j60 = np.argmin(np.abs( Cv_60 - np.max(Cv_60)  ))
j80 = np.argmin(np.abs( Cv_80 - np.max(Cv_80)  ))
j100 = np.argmin(np.abs( Cv_100 - np.max(Cv_100)  ))

Tc_40 = T40[j40]
Tc_60 = T60[j60]
Tc_80 = T80[j80]
Tc_100 = T100[j100]

print(" L :    |    40     |    60    |    80    |    100    |")
print("Tc(Cv) : ", Tc_40, Tc_60, Tc_80, Tc_100 )






#     ---------- PROBLEM 9 ----------



Tc_list = np.array([Tc_40, Tc_60, Tc_80, Tc_100])


# function fitted to data to find Tc(L=inf)
def f(x,a,b):
    return a/x + b

a,b = curve_fit(f,L,Tc_list)  #fitting


print("a and Tc(L=inf) :", a[0],a[1])

# plotting the fitted curve with the data
L_range = np.linspace(2,150,1000)

plt.plot(L, Tc_list,"o", label="data")
plt.plot(L_range, f(L_range,a[0],a[1]), label=f"Tc(L): a={a[0]:.4f}, " + r"$L_{\infty}$" + f"={a[1]:.4f}")
plt.ylabel(r"$T_c$ $[J/k_B]$")
plt.xlabel(r"Lattice size, L")
plt.legend()
plt.savefig("figs/p8_Tc_fit.pdf")
plt.show()