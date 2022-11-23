import numpy as np
import matplotlib.pyplot as plt


#     ---------- PROBLEM 4 ----------



N = 2*2  # number of spins

#   --- analytical ---
T = 1                  # temperature
beta = 1/(T)   



E2x2 = np.array([-8, 0, 0, 0, -8, 8])    # array of energies
e2x2 = E2x2/N
D2x2 = np.array([ 1, 4, 4, 4, 1, 2])   # array of degenerasi

M2x2 = np.array([4, 2, 0, -2, -4, 0])  # magentization
m2x2= M2x2/N
Z2x2 = 0
for i in range(len(E2x2)): 
    Z2x2 += np.exp(-beta*E2x2[i])*D2x2[i]   # finding the partition function


P2x2 = 1/Z2x2 * np.exp(-beta*E2x2)*D2x2  #percent chance for each state

avrg_E_a = np.sum(E2x2*P2x2)      # averag E
avrg_E2_a = np.sum(E2x2**2*P2x2)  # average E^2

avrg_M_a = np.sum(np.abs(M2x2)*P2x2) #average M
avrg_M2_a = np.sum(M2x2**2*P2x2)  #average M^2


print("<e>_analytical=", avrg_E_a/N)
print("<e^2>_analytical=", avrg_E2_a/N)
print()
print("<|m|>_analytical=", avrg_M_a/N)
print( "<m^2>_analytical=", avrg_M2_a/N)

Cv_a = 1/(N*T**2)*(avrg_E2_a - avrg_E_a**2)  # C_V and chi
Xi_a = 1/(N*T)*(avrg_M2_a - avrg_M_a**2)

print()
print("Cv_a=", Cv_a)
print( "chi_a=", Xi_a) 
print()


#   --- numerical ---


EM = np.loadtxt("2x2_n=10^6_all.txt")


start = 0
end = 10**6
# print(EM[0,1])
E = EM[start:end+1,0]
M = EM[start:end+1,1]

NN = 6                    # calculating the values at 10^(j+1) MC cycles
avrg_E_l, avrg_M_l, avrg_E2_l, avrg_M2_l, Cv_l, Xi_l = np.zeros(NN), np.zeros(NN), np.zeros(NN), np.zeros(NN), np.zeros(NN), np.zeros(NN)
for j in range(NN):
    i = 10**(j+1)
    avrg_E_l[j] = np.mean(E[:i])
    avrg_M_l[j] = np.mean(np.abs(M[:i]))
    avrg_E2_l[j] = np.mean((E[:i])**2)
    avrg_M2_l[j] = np.mean((M[:i])**2)
    Cv_l[j] = 1/(N*T**2)*(avrg_E2_l[j] - avrg_E_l[j]**2)
    Xi_l[j] = 1/(N*T)*(avrg_M2_l[j] - avrg_M_l[j]**2)


print( f"  n    |   <e>    |  <e^2>  |   <m> | <m^2>  | Cv   | Xi    |"  )
for i in range(NN):
    print(f"{10**(i+1):5.5g}  |" ,f"{avrg_E_l[i]/N:6.6f}",  f"{avrg_E2_l[i]/N**2:6.6f}", f"{avrg_M_l[i]/N:6.6f}",  f"{avrg_M2_l[i]/N**2:6.6f}", f"{Cv_l[i]:6.6f}",  f"{Xi_l[i]:6.6f}"  )