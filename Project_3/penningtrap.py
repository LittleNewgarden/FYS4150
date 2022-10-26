from cProfile import label
from unittest.util import three_way_cmp
import numpy as np
import matplotlib.pyplot as plt



# finding wz, max R and min R for the analytical solution

q = 1
V0 = 2.41*10**6
B0 = 9.65*10
m =  40.078
z0 = 20
d = 500
wz = np.sqrt(2*q*V0/(m*d**2))
w0 = q*B0/m
w_p = (w0+np.sqrt(w0**2 -2*wz**2))/2
w_m = (w0-np.sqrt(w0**2 -2*wz**2))/2
v0 = 25
x0 = 20
A_p = (v0+w_m*x0)/(w_m - w_p)
A_m = -(v0+w_p*x0)/(w_m - w_p)

print("max R:", A_p+A_m)
print("min R:", np.abs(A_p-A_m))

print("Wz analytical:",wz)



# 8 point one


R_x,R_y,R_z = np.loadtxt("r_n_iterations=4000_method=RK4_number_particles=1_particle=0_interactions=1.txt",unpack=True)

t = np.linspace(0,50,4001)

# plotting the z as function of time
plt.plot(t,R_z)
plt.xlabel(r"time [$\mu s]$")
plt.ylabel(r"z [$\mu m$]")
plt.savefig("t_z.pdf")
plt.show()


# printing wz found numerically
j1 = np.argmin(np.abs(t-5))
j2 = np.argmin(np.abs(t-15))
j3 = np.argmin(np.abs(t-25))

p1 = np.where(R_z==np.max(R_z[j1:j2]))
p2 = np.where(R_z==np.max(R_z[j2:j3]))

print("w_z numerical:", 2*np.pi/(t[p2]-t[p1]))


# finding r min and r max numerically
abs_r = np.sqrt(R_x**2 + R_y**2)
r_min = np.min(abs_r)
r_max = np.max(abs_r)
print("r_min numerical:",r_min)
print("r_max numerical:",r_max)



# point two

# with interactions xy-plot
R_x1_1, R_y1_1, R_z1_1 = np.loadtxt("r_n_iterations=4000_method=RK4_number_particles=2_particle=0_interactions=1.txt",unpack=True)
R_x2_1, R_y2_1, R_z2_1 = np.loadtxt("r_n_iterations=4000_method=RK4_number_particles=2_particle=1_interactions=1.txt",unpack=True)


plt.plot(R_x1_1,R_y1_1,label="Particle 1")
plt.plot(R_x2_1,R_y2_1,label="Particle 2")
plt.scatter(R_x1_1[0],R_y1_1[0])
plt.scatter(R_x2_1[0],R_y2_1[0])
plt.xlabel(r"x [$\mu m$]")
plt.ylabel(r"y [$\mu m$]")
plt.axis('equal')
plt.legend()
plt.savefig("x_y_interactions.pdf")
plt.show()


# Without interactions xy-plot
R_x1_0, R_y1_0, R_z1_0 = np.loadtxt("r_n_iterations=4000_method=RK4_number_particles=2_particle=0_interactions=0.txt",unpack=True)
R_x2_0, R_y2_0, R_z2_0 = np.loadtxt("r_n_iterations=4000_method=RK4_number_particles=2_particle=1_interactions=0.txt",unpack=True)


plt.plot(R_x1_0,R_y1_0,label="Particle 1")
plt.plot(R_x2_0,R_y2_0,label="Particle 2")
plt.scatter(R_x1_0[0],R_y1_0[0])
plt.scatter(R_x2_0[0],R_y2_0[0])
plt.xlabel(r"x [$\mu m$]")
plt.ylabel(r"y [$\mu m$]")
plt.axis('equal')
plt.legend()
plt.savefig("x_y_no_interactions.pdf")
plt.show()




# point three

V_x1_1, V_y1_1, V_z1_1 = np.loadtxt("v_n_iterations=4000_method=RK4_number_particles=2_particle=0_interactions=1.txt",unpack=True)
V_x2_1, V_y2_1, V_z2_1 = np.loadtxt("v_n_iterations=4000_method=RK4_number_particles=2_particle=1_interactions=1.txt",unpack=True)

V_x1_0, V_y1_0, V_z1_0 = np.loadtxt("v_n_iterations=4000_method=RK4_number_particles=2_particle=0_interactions=0.txt",unpack=True)
V_x2_0, V_y2_0, V_z2_0 = np.loadtxt("v_n_iterations=4000_method=RK4_number_particles=2_particle=1_interactions=0.txt",unpack=True)


# with interactions x,vx plot and z,vz plot
plt.plot(R_x1_1,V_x1_1,label="Particle 1")
plt.plot(R_x2_1,V_x2_1,label="Particle 2")
plt.scatter(R_x1_1[0],V_x1_1[0])
plt.scatter(R_x2_1[0],V_x2_1[0])
plt.xlabel(r"$x$ [$\mu m$]")
plt.ylabel(r"$v_x$ [$\mu m/\mu s$]")
plt.axis('equal')
plt.legend()
plt.savefig("x_v_interactions.pdf")
plt.show()

plt.plot(R_z1_1,V_z1_1,label="Particle 1")
plt.plot(R_z2_1,V_z2_1,label="Particle 2")
plt.scatter(R_z1_1[0],V_z1_1[0])
plt.scatter(R_z2_1[0],V_z2_1[0])
plt.xlabel(r"$z$ [$\mu m$]")
plt.ylabel(r"$v_z$ [$\mu m/\mu s$]")
plt.axis('equal')
plt.legend()
plt.savefig("z_v_interactions.pdf")
plt.show()


# without interactosx,vx plot and z,vz plot
plt.plot(R_x1_0,V_x1_0,label="Particle 1")
plt.plot(R_x2_0,V_x2_0,label="Particle 2")
plt.scatter(R_x1_0[0],V_x1_0[0])
plt.scatter(R_x2_0[0],V_x2_0[0])
plt.xlabel(r"$x$ [$\mu m$]")
plt.ylabel(r"$v_x$ [$\mu m/\mu s$]")
plt.axis('equal')
plt.legend()
plt.savefig("x_v_no_interactions.pdf")
plt.show()

plt.plot(R_z1_0,V_z1_0,label="Particle 1")
plt.plot(R_z2_0,V_z2_0,label="Particle 2")
plt.scatter(R_z1_0[0],V_z1_0[0])
plt.scatter(R_z2_0[0],V_z2_0[0])
plt.xlabel(r"$z$ [$\mu m$]")
plt.ylabel(r"$v_z$ [$\mu m/\mu s$]")
plt.axis('equal')
plt.legend()
plt.savefig("z_v_no_interactions.pdf")
plt.show()





# point four, 3d plot with and without interactions
ax = plt.axes(projection = "3d")
ax.plot(R_x1_1, R_y1_1, R_z1_1,label="Particle 1")
ax.plot(R_x2_1, R_y2_1, R_z2_1,label="Particle 2")
ax.scatter(R_x1_1[0], R_y1_1[0], R_z1_1[0])
ax.scatter(R_x2_1[0], R_y2_1[0], R_z2_1[0])
ax.set_xlabel(r"$x$ [$\mu m$]")
ax.set_ylabel(r"$y$ [$\mu m$]")
ax.set_zlabel(r"$z$ [$\mu m$]")
plt.legend()
plt.savefig("3d_interactions.pdf")
plt.show()


ax = plt.axes(projection = "3d")
ax.plot(R_x1_0, R_y1_0, R_z1_0,label="Particle 1")
ax.plot(R_x2_0, R_y2_0, R_z2_0,label="Particle 2")
ax.scatter(R_x1_0[0], R_y1_0[0], R_z1_0[0])
ax.scatter(R_x2_0[0], R_y2_0[0], R_z2_0[0])
ax.set_xlabel(r"$x$ [$\mu m$]")
ax.set_ylabel(r"$y$ [$\mu m$]")
ax.set_zlabel(r"$z$ [$\mu m$]")
plt.legend()
plt.savefig("3d_no_interactions.pdf")
plt.show()





# point 5


n_l = np.array([4000,8000,16000,32000])

R_Euler = []
R_RK4 = []


for i in range(len(n_l)):
    R_Euler.append(np.loadtxt(f"r_n_iterations={n_l[i]}_method=Euler_number_particles=1_particle=0_interactions=1.txt",unpack=False))
    R_RK4.append(np.loadtxt(f"r_n_iterations={n_l[i]}_method=RK4_number_particles=1_particle=0_interactions=1.txt",unpack=False))



# analytical solution x,,y,z
def R_analytical(t):
    q = 1
    V0 = 2.41*10**6
    B0 = 9.65*10
    m =  40.078
    z0 = 20
    d = 500
    wz = np.sqrt(2*q*V0/(m*d**2))
    w0 = q*B0/m
    w_p = (w0+np.sqrt(w0**2 -2*wz**2))/2
    w_m = (w0-np.sqrt(w0**2 -2*wz**2))/2
    v0 = 25
    x0 = 20
    A_p = (v0+w_m*x0)/(w_m - w_p)
    A_m = -(v0+w_p*x0)/(w_m - w_p)

    x = A_p*np.cos(w_p*t) + A_m*np.cos(w_m*t)
    y = -A_p*np.sin(w_p*t) - A_m*np.sin(w_m*t)
    z = z0*np.cos(wz*t)
    return np.array([x, y, z])


#different time axis
t1 = np.linspace(0,50,4001)
t2 = np.linspace(0,50,8001)
t3 = np.linspace(0,50,16001)
t4 = np.linspace(0,50,32001)

t_l = [t1,t2,t3,t4]


error = np.full((4,32001), np.nan)
error_E = np.full((4,32001), np.nan)


delta_RK4 = np.full((4,32001), np.nan)
delta_Euler = np.full((4,32001), np.nan)

#findinf delta and absolute error for rk4 and forward euler
for i in range(len(n_l)):
    for j in range(len(t_l[i])):
        delta_RK4[i,j] = np.linalg.norm(  R_analytical(t_l[i][j]) - R_RK4[i][j,:]  )
        delta_Euler[i,j] = np.linalg.norm(  R_analytical(t_l[i][j]) - R_Euler[i][j,:]  )

        error[i,j] = delta_RK4[i,j]/np.linalg.norm(R_analytical(t_l[i][j]))
        error_E[i,j] = delta_Euler[i,j]/np.linalg.norm(R_analytical(t_l[i][j]))
        
        

# plotting absolute error of rk4 and forward euler
for i in range(len(t_l)):
    b = np.isnan(error[i])==False
    plt.plot(t_l[i], error[i,b],label=f"n = {len(t_l[i])-1}")
plt.yscale('log')
plt.xlabel(r"time [$\mu s$]")
plt.ylabel(r"relative error")
plt.legend()
plt.savefig("RK4_rel_err.pdf")
plt.show()

for i in range(len(t_l)):
    b = np.isnan(error[i])==False
    plt.plot(t_l[i], error_E[i,b],label=f"n = {len(t_l[i])-1}")
plt.yscale('log')
plt.xlabel(r"time [$\mu s$]")
plt.ylabel(r"relative error")
plt.legend()
plt.savefig("Euler_rel_err.pdf")
plt.show()





# point 6 calculatin rhe error convergence rate


d_max_RK4 = np.zeros(4)
d_max_Euler = np.zeros(4)
h = np.zeros(4)

for i in range(4):
    d_max_RK4[i] = np.nanmax(  delta_RK4[i]  )
    d_max_Euler[i] = np.nanmax(  delta_Euler[i]  )
    h[i] = t_l[i][1] - t_l[i][0]


r_err_k_RK4 = np.log(d_max_RK4[1:]/d_max_RK4[:-1])/np.log(h[1:]/h[:-1])
r_err_RK4 = 1/3*np.sum(r_err_k_RK4)

r_err_k_Euler = np.log(d_max_Euler[1:]/d_max_Euler[:-1])/np.log(h[1:]/h[:-1])
r_err_Euler = 1/3*np.sum(r_err_k_Euler)

 
print(" error convergence rate RK4", r_err_RK4)
print(" error convergence rate Euler", r_err_Euler)



##########  9 1  ############### plotting the fraction of particles left with no interaction



N_end = np.loadtxt(f"100_particles_500ms_interactions=0.txt")
N_tot = 100
steps =  int((2.5-0.2)/0.02) + 1
print("steps",steps)
Wv = np.linspace(0.2,2.5,steps)
f_l = [0.1,0.4,0.7] #different amplitudes
for i in range(3):
    plt.plot(Wv,N_end[i,:]/N_tot, label=f"f = {f_l[i]}")
plt.xlabel(r"$\omega_V$ [$1/\mu s$]")
plt.ylabel(r"N/N_{tot}")
plt.legend()
plt.savefig("Wv_N_full.pdf")
plt.show()




############   9  2   ############ plotting the fraction of particles left with and without interactions for limited range


N_end = np.loadtxt(f"100_particles_500ms_interactions=0_zoomed_inn_n=50_a.txt")
N_end_inter = np.loadtxt(f"100_particles_500ms_interactions=1_zoomed_inn_n=50_a.txt")

N_end2 = np.loadtxt(f"100_particles_500ms_interactions=0_zoomed_inn_n=50_b.txt")
N_end_inter2 = np.loadtxt(f"100_particles_500ms_interactions=1_zoomed_inn_n=50_b.txt")

N_no_interactions = (N_end + N_end2)/2

N_interactions = (N_end_inter + N_end_inter2)/2


steps = 50
Wv = np.linspace(2.12,2.22,steps)


plt.plot(Wv, N_no_interactions/N_tot, label="No interactions")
plt.plot(Wv, N_interactions/N_tot,label="Interactions")
plt.xlabel(r"$\omega_V$ [$1/\mu s$]")
plt.ylabel(r"$N/N_{tot}$")
plt.legend()
plt.savefig("Wv_N.pdf")
plt.show()

