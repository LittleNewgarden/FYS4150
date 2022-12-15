import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support

quantity_support()
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(10, 8)) # Larger figure sizes
plt.rc('font', size=18)



# ----------------  P7  -------------

# number of points, M. number of timesteps, Tn. end time, T 

M = 201

Tn = 321
Tend = 0.008

T = np.linspace(0,Tend, Tn)


# V = 0 at all points
P_1 = np.loadtxt("task7_1.txt")

# double slit
P_2 = np.loadtxt("task7_2.txt")

p_1 = np.zeros((Tn, M, M))
p_2 = np.zeros((Tn, M, M))

# putting the cubes in a more easely eccesable format
for i in range(Tn):
    p_1[i,:,:] = P_1[M*i:M*(i+1), :]
    p_2[i,:,:] = P_2[M*i:M*(i+1), :]


# plotting the error for no wall and double slit

p1_tot = np.sum(np.sum(p_1, axis=1 ), axis=1)
plt.plot(T,p1_tot - 1)
plt.ylabel(r"$\sum_{i,j}(p_{i,j}) - 1$")
plt.xlabel("t")
plt.savefig("figures/p7_1.pdf")
plt.show()


p2_tot = np.sum(np.sum(p_2, axis=1 ), axis=1)
plt.plot(T,p2_tot - 1)
plt.ylabel(r"$\sum_{i,j}(p_{i,j}) - 1$")
plt.xlabel("t")
plt.savefig("figures/p7_2.pdf")
plt.show()











#  --------------------   P8   -------------------

# number of points, M. number of timesteps, Tn. end time, T 

M = 201

Tn = 81
Tend = 0.002

T = np.linspace(0,Tend, Tn)


P = np.loadtxt("task8.txt")
U_real = np.loadtxt("task8_real.txt")
U_imag = np.loadtxt("task8_imag.txt")

p = np.zeros((Tn, M, M))
u_real = np.zeros((Tn, M, M))
u_imag = np.zeros((Tn, M, M))

# putting the cubes in a more easely eccesable format
for i in range(Tn):
    p[i,:,:] = P[M*i:M*(i+1), :]
    u_real[i,:,:] = U_real[M*i:M*(i+1), :]
    u_imag[i,:,:] = U_imag[M*i:M*(i+1), :]


# finding the index of t = 0.001
j = np.argmin( np.abs(T-0.001))

indx = [0, j, -1]

# cahnging the x and y axis in imshow
x = np.linspace(0,1,5)
grid = np.linspace(0,200,5)



# plotting P for all three times

for i in indx:
    im = plt.imshow(p[i], origin='lower') 
    plt.colorbar(im, label = f"p(t = {T[i]})")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks(grid, x)
    plt.yticks(grid, x)
    plt.savefig(f"figures/p8_p_{i}.pdf")
    plt.show()

# plotting real(U) for all three times


for i in indx:
    im = plt.imshow(u_real[i], origin='lower')
    plt.colorbar(im, label = f"Re(u(t = {T[i]}))")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks(grid, x)
    plt.yticks(grid, x)
    plt.savefig(f"figures/p8_u_real_{i}.pdf")
    plt.show()


# plotting imag(U) for all three times

for i in indx:
    im = plt.imshow(u_imag[i], origin='lower')
    plt.colorbar(im, label = f"Im(u(t = {T[i]}))")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks(grid, x)
    plt.yticks(grid, x)
    plt.savefig(f"figures/p8_u_imag_{i}.pdf")
    plt.show()














#   -------     P9   -------


# number of points, M. number of timesteps, Tn. end time, T 


M = 201

Tn = 81
Tend = 0.002

T = np.linspace(0,Tend, Tn)

# single slit
P_s = np.loadtxt("task9_single.txt")

# double slit
# P_d = P #np.loadtxt("task9_double.txt")

# triple slit
P_t = np.loadtxt("task9_triple.txt")

# single slit
p_s = np.zeros((Tn, M, M))

# double slit
p_d = p  #np.zeros((Tn, M, M))

# triple slit
p_t = np.zeros((Tn, M, M))

# putting the cubes in a more easely eccesable format
for i in range(Tn):
    p_s[i,:,:] = P_s[M*i:M*(i+1), :]
    # p_d[i,:,:] = P[M*i:M*(i+1), :]
    p_t[i,:,:] = P_t[M*i:M*(i+1), :]






# finding the index of the screen 
screen_x = 0.8
x2 = np.linspace(0,1,M)

j_x = np.argmin(np.abs(x2 - screen_x))



# creating a line accros the imshow to visualise the screen 
screen_x = np.ones(M)*j_x
screen_y = np.linspace(0,M-1, M)


# plotting single slit P for t = 0.002

im = plt.imshow(p_s[-1], origin='lower') #, cmap=plt.get_cmap('gist_gray'))
plt.colorbar(im, label = f"p(t = {T[-1]})")
plt.plot(screen_x, screen_y, color="black", label="screen")
plt.xlabel("x")
plt.ylabel("y")
plt.xticks(grid, x)
plt.yticks(grid, x)
plt.legend()
plt.savefig(f"figures/p9_im_1.pdf")
plt.show()


# plotting normalised single slit p along y axis for x = 0.8

plt.plot(x2, p_s[-1, :, j_x]/np.sum(p_s[-1, :, j_x]))
plt.xlabel("y")
plt.ylabel("p(y|x = 0.8; t=0.002)")
plt.savefig("figures/p9_1.pdf")
plt.show()


# plotting double slit P for t = 0.002

im = plt.imshow(p_d[-1], origin='lower') #, cmap=plt.get_cmap('gist_gray'))
plt.colorbar(im, label = f"p(t = {T[-1]})")
plt.plot(screen_x, screen_y, color="black", label="screen")
plt.xlabel("x")
plt.ylabel("y")
plt.xticks(grid, x)
plt.yticks(grid, x)
plt.legend()
plt.savefig(f"figures/p9_im_2.pdf")
plt.show()


# plotting normalised double slit p along y axis for x = 0.8

plt.plot(x2, p_d[-1, :, j_x]/np.sum(p_d[-1, :, j_x]))
plt.xlabel("y")
plt.ylabel("p(y|x = 0.8; t=0.002)")
plt.savefig("figures/p9_2.pdf")
plt.show()


# plotting triple slit P for t = 0.002

im = plt.imshow(p_t[-1], origin='lower') #, cmap=plt.get_cmap('gist_gray'))
plt.plot(screen_x, screen_y, color="black", label="screen")
plt.colorbar(im, label = f"p(t = {T[-1]})")
plt.xlabel("x")
plt.ylabel("y")
plt.xticks(grid, x)
plt.yticks(grid, x)
plt.legend()
plt.savefig(f"figures/p9_im_3.pdf")
plt.show()


# plotting normalised triple slit p along y axis for x = 0.8

plt.plot(x2, p_t[-1, :, j_x]/np.sum(p_t[-1, :, j_x]))
plt.xlabel("y")
plt.ylabel("p(y|x = 0.8; t=0.002)")
plt.savefig("figures/p9_3.pdf")
plt.show()


