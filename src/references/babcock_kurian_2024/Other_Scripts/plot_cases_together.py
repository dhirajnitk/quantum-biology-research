from numpy import array, complex64
import numpy as np
from matplotlib import pyplot as plt
import os

h = 6.626e-27                # cm^2 g s^-1
c = 2.998e10                 # cm s^-1
e0 = 35716.65                # cm^-1        (energy)
k0 = 2.24e-3                 # A^-1         (inverse length)
musq = 181224                # A^3 cm^-1 
gamma = 2.73e-3              # cm^-1        (energy)
w0 = c*k0*(1e8)              # s^-1         (angular frequency)
ns = 104

AX_LAB_SIZE = 25
AX_TICK_SIZE = 20
AX_LEGEND_SIZE = 15

def index_of_SRS(eigenvalues):
	Energ_eigenvalues = np.real(eigenvalues)
	Gamma_eigenvalues = np.imag(eigenvalues) / (-gamma/2)
	return np.where(np.sort(Energ_eigenvalues) == Energ_eigenvalues[np.where(Gamma_eigenvalues == np.max(Gamma_eigenvalues))[0][0]])[0][0]


with open("Eigenvalues/realistic_MT_100_spiral_approximate_eigvals.txt", 'r') as f:
	real_approx_eigvals_str = f.readlines()
	real_approx_eigvals = eval(" ".join(real_approx_eigvals_str).replace("\n", ""))


with open("Eigenvalues/realistic_MT_100_spiral_numerical_eigvals.txt", 'r') as f:
	real_num_eigvals_str = f.readlines()
	real_num_eigvals = eval(" ".join(real_num_eigvals_str).replace("\n", ""))


with open("Eigenvalues/z_proj_realistic_MT_100_spiral_approximate_eigvals.txt", 'r') as f:
	z_proj_approx_eigvals_str = f.readlines()
	z_proj_approx_eigvals = eval(" ".join(z_proj_approx_eigvals_str).replace("\n", ""))


with open("Eigenvalues/z_proj_realistic_MT_100_spiral_numerical_eigvals.txt", 'r') as f:
	z_proj_num_eigvals_str = f.readlines()
	z_proj_num_eigvals = eval(" ".join(z_proj_num_eigvals_str).replace("\n", ""))

with open("Eigenvalues/perfect_cylinder_eigenvalues_100_spiral.txt", 'r') as f:
	perf_cylinder_eigvals_str = f.readlines()
	perf_cylinder_eigvals = eval(" ".join(perf_cylinder_eigvals_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_approx_str = f.readlines()
	perf_ring_eigvals_approx = eval(" ".join(perf_ring_eigvals_approx_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_numerical_diag_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_numerical_str = f.readlines()
	perf_ring_eigvals_numerical = eval(" ".join(perf_ring_eigvals_numerical_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_without_LR_int_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_approx_without_LR_str = f.readlines()
	perf_ring_eigvals_approx_without_LR = eval(" ".join(perf_ring_eigvals_approx_without_LR_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_numerical_diag_without_LR_int_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_numerical_without_LR_str = f.readlines()
	perf_ring_eigvals_numerical_without_LR = eval(" ".join(perf_ring_eigvals_numerical_without_LR_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_noise_5_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_noise_5_str = f.readlines()
	perf_ring_eigvals_noise_5 = eval(" ".join(perf_ring_eigvals_noise_5_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_noise_10_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_noise_10_str = f.readlines()
	perf_ring_eigvals_noise_10 = eval(" ".join(perf_ring_eigvals_noise_10_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_noise_20_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_noise_20_str = f.readlines()
	perf_ring_eigvals_noise_20 = eval(" ".join(perf_ring_eigvals_noise_20_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_noise_50_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_noise_50_str = f.readlines()
	perf_ring_eigvals_noise_50 = eval(" ".join(perf_ring_eigvals_noise_50_str).replace("\n", ""))

with open("Eigenvalues/perfect_rings_PD_eigenvalues_approx_diag_noise_100_100_spiral.txt", 'r') as f:
	perf_ring_eigvals_noise_100_str = f.readlines()
	perf_ring_eigvals_noise_100 = eval(" ".join(perf_ring_eigvals_noise_100_str).replace("\n", ""))

#print(np.sort(np.imag(real_approx_eigvals) / (-gamma/2))[-1])
"""
op = 0.8

plt.figure(figsize=(12,8))
plt.rc('axes', labelsize=AX_LAB_SIZE) 
plt.rc('xtick', labelsize=AX_TICK_SIZE)
plt.rc('ytick', labelsize=AX_TICK_SIZE) 
plt.rc('legend', fontsize=AX_LEGEND_SIZE) 
plt.rc('axes', titlesize=AX_LAB_SIZE)
pcsp = plt.scatter(np.real(perf_cylinder_eigvals), np.imag(perf_cylinder_eigvals) / (-gamma/2), c='g', alpha=op, label="Diag. of MT structure as perfect rings with tensor product basis")
rnsp = plt.scatter(np.real(real_num_eigvals), np.imag(real_num_eigvals) / (-gamma/2), c='b', alpha=op, label="Exact diag. of realistic MT dipoles (model 1)")
rasp = plt.scatter(np.real(real_approx_eigvals), np.imag(real_approx_eigvals) / (-gamma/2), c='c', alpha=op, label="Approximate diag. of realistic MT dipoles (model 2)")
znsp = plt.scatter(np.real(z_proj_num_eigvals), np.imag(z_proj_num_eigvals) / (-gamma/2), c='r', alpha=op, label="Exact diag. of z-projected MT dipoles (model 1)")
zasp = plt.scatter(np.real(z_proj_approx_eigvals), np.imag(z_proj_approx_eigvals) / (-gamma/2), c='m', alpha=op, label="Approximate diag. of z-projected MT dipoles (model 2)")

plt.xlim(-150, 150)
plt.ylim(-10, 600)
plt.xlabel("Energy [cm^-1]")
plt.ylabel("Γ / γ")
plt.legend(handles=[rnsp, rasp, znsp, zasp, pcsp])

plt.figure(figsize=(12,8))

#finite_continuous_cylinder_energies = np.array([-5.83608,-5.66337,-5.61505,-5.33854,-5.12015,-4.74752,-4.44607,-4.08085,-3.72278,-3.46218,-3.07593,-2.9826,-2.68071,-2.59495,-2.54041,-2.50538,-2.50405,-2.47583,-2.38952,-2.31488,-2.31421,-2.28982,-2.28262,-2.24879,-2.227,-2.21215,-2.20859,-2.09898,-2.08425,-1.99871,-1.97155,-1.97056,-1.94442,-1.93727,-1.93669,-1.92371,-1.85501,-1.84218,-1.83868,-1.82653,-1.81292,-1.80169,-1.75028,-1.63914,-1.51857,-1.41553,-1.34784,-1.31816,-1.31762,-1.31467,-1.30843,-1.27763,-1.2281,-1.17289,-1.16699,-1.15873,-1.149,-1.135,-1.12861,-1.11266,-1.10768,-1.10412,-1.04177,-0.977218,-0.925274,-0.894593,-0.891297,-0.887667,-0.885754,-0.884811,-0.867298,-0.850025,-0.84645,-0.839687,-0.833869,-0.828555,-0.812068,-0.806931,-0.79983,-0.794471,-0.793022,-0.756125,-0.71389,-0.682714,-0.672386,-0.66951,-0.668563,-0.668111,-0.667604,-0.667054,-0.664947,-0.655855,-0.654146,-0.652075,-0.637061,-0.63616,-0.623618,-0.62103,-0.620819,-0.614802,-0.590136])

plt.rc('axes', labelsize=AX_LAB_SIZE) 
plt.rc('xtick', labelsize=AX_TICK_SIZE)
plt.rc('ytick', labelsize=AX_TICK_SIZE) 
plt.rc('legend', fontsize=AX_LEGEND_SIZE) 
plt.rc('axes', titlesize=AX_LAB_SIZE)
rnp = plt.plot(np.array(range(10400)), np.sort(np.real(real_num_eigvals)), c='b', label="Exact diag. of realistic MT dipoles (model 1)")
rap = plt.plot(np.array(range(10400)), np.sort(np.real(real_approx_eigvals)), c='c', label="Approximate diag. of realistic MT dipoles (model 2)")
znp = plt.plot(np.array(range(10400)), np.sort(np.real(z_proj_num_eigvals)), c='r', label="Exact diag. of z-projected MT dipoles (model 1)")
zap = plt.plot(np.array(range(10400)), np.sort(np.real(z_proj_approx_eigvals)), c='m', label="Approximate diag. of z-projected MT dipoles (model 2)")
pcp = plt.plot(np.array(range(10400)), np.sort(np.real(perf_cylinder_eigvals)), c='g', label="Diag. of MT structure as perfect rings with tensor product basis")
#fcp = plt.plot(np.array(range(np.size(finite_continuous_cylinder_energies))), finite_continuous_cylinder_energies, c='g', label="Finite cylinder with continuous distribution of z-oriented dipoles")

SRS_indices = np.array([index_of_SRS(real_num_eigvals), index_of_SRS(real_approx_eigvals), index_of_SRS(z_proj_num_eigvals), index_of_SRS(z_proj_approx_eigvals), index_of_SRS(perf_cylinder_eigvals)])
max_SRS = np.array([np.sort(np.real(real_num_eigvals))[SRS_indices[0]], np.sort(np.real(real_approx_eigvals))[SRS_indices[1]], np.sort(np.real(z_proj_num_eigvals))[SRS_indices[2]],np.sort(np.real(z_proj_approx_eigvals))[SRS_indices[3]], np.sort(np.real(perf_cylinder_eigvals))[SRS_indices[4]]])

plt.scatter(SRS_indices, max_SRS, c='k', label="Maximally superradiant states")
plt.ylim(-150, 150)
plt.xlabel("State index (sorted)")
plt.ylabel("Energy [cm^-1]")
plt.legend()

plt.figure(figsize=(12,8))

#plt.scatter(np.real(perf_ring_eigvals_approx_without_LR), np.imag(perf_ring_eigvals_approx_without_LR) / (-gamma/2), c='r')
#plt.scatter(np.real(perf_ring_eigvals_numerical_without_LR), np.imag(perf_ring_eigvals_numerical_without_LR) / (-gamma/2), c='g')
#plt.scatter(np.real(perf_ring_eigvals_numerical), np.imag(perf_ring_eigvals_numerical) / (-gamma/2), c='g')

#imdiff = np.sort(np.imag(perf_ring_eigvals_approx) / (-gamma/2)) - np.sort(np.imag(perf_ring_eigvals_numerical) / (-gamma/2))
#rediff = np.sort(np.real(perf_ring_eigvals_approx)) - np.sort(np.real(perf_ring_eigvals_numerical))

#plt.scatter(np.array(range(10400)), imdiff, c='red', s=3)
#plt.scatter(np.array(range(10400)), rediff, c='black', s=3)

# Approx. vs. numerical without LR interactions
plt.figure(figsize=(12,8))

plt.scatter(np.real(perf_ring_eigvals_approx_without_LR), np.imag(perf_ring_eigvals_approx_without_LR) / (-gamma/2), c='red')
plt.scatter(np.real(perf_ring_eigvals_numerical_without_LR), np.imag(perf_ring_eigvals_numerical_without_LR) / (-gamma/2), c='black')

# Eigenvalues with differing amounts of thermal noise.
plt.figure(figsize=(12,8))

n100 = plt.scatter(np.real(perf_ring_eigvals_noise_100), np.imag(perf_ring_eigvals_noise_100) / (-gamma/2), c='blue', label="W = 100")
n50 = plt.scatter(np.real(perf_ring_eigvals_noise_50), np.imag(perf_ring_eigvals_noise_50) / (-gamma/2), c='green', label="W = 50")
n20 = plt.scatter(np.real(perf_ring_eigvals_noise_20), np.imag(perf_ring_eigvals_noise_20) / (-gamma/2), c='yellow', label="W = 20")
n10 = plt.scatter(np.real(perf_ring_eigvals_noise_10), np.imag(perf_ring_eigvals_noise_10) / (-gamma/2), c='orange', label="W = 10")
n5 = plt.scatter(np.real(perf_ring_eigvals_noise_5), np.imag(perf_ring_eigvals_noise_5) / (-gamma/2), c='red', label="W = 5")

plt.legend(handles=[n5,n10,n20,n50,n100])

# Approx. eigenvalues (for PD/z-oriented case) vs numerical eigenvals with LR interactions
plt.figure(figsize=(12,8))

plt.scatter(np.real(perf_ring_eigvals_approx), np.imag(perf_ring_eigvals_approx) / (-gamma/2), c='red')
plt.scatter(np.real(perf_ring_eigvals_numerical), np.imag(perf_ring_eigvals_numerical) / (-gamma/2), c='black')

# Residuals of approx. eigvals vs numerical eigvals
plt.figure(figsize=(12,8))

imdiff = np.sort(np.imag(perf_ring_eigvals_approx) / (-gamma/2)) - np.sort(np.imag(perf_ring_eigvals_numerical) / (-gamma/2))
rediff = np.sort(np.real(perf_ring_eigvals_approx)) - np.sort(np.real(perf_ring_eigvals_numerical))
plt.scatter(np.array(range(10400)), imdiff, c='red', s=3)
plt.scatter(np.array(range(10400)), rediff, c='black', s=3)
"""
# Reproduction of fig. S5 in preprint (how max. superradiance varies with thermal noise)

w_to_ind = {
	0: 0,
	5: 1,
	10: 2,
	20: 3,
	50: 4,
	100: 5
}

directory = "Eigenvalues/"
seg_10 = np.zeros((6))
seg_20 = np.zeros((6))
seg_50 = np.zeros((6))
seg_80 = np.zeros((6))

for filename in os.listdir(directory):
	f = os.path.join(directory, filename)
	if os.path.isfile(f):
		if f[-4:] == ".npy" and "eigenvalues" in f:
			fsplit = f.split("_")
			if fsplit[3] == "10":
				seg_10[w_to_ind[int(fsplit[6].split(".")[0])]] = np.max(np.imag(np.load(f)) / (-gamma/2))
			
			elif fsplit[3] == "20":
				seg_20[w_to_ind[int(fsplit[6].split(".")[0])]] = np.max(np.imag(np.load(f)) / (-gamma/2))
			
			elif fsplit[3] == "50":
				seg_50[w_to_ind[int(fsplit[6].split(".")[0])]] = np.max(np.imag(np.load(f)) / (-gamma/2))

			elif fsplit[3] == "80":
				seg_80[w_to_ind[int(fsplit[6].split(".")[0])]] = np.max(np.imag(np.load(f)) / (-gamma/2))

plt.figure(figsize=(12,8))

max_decay = np.max(np.imag(perf_ring_eigvals_numerical) / (-gamma/2))

x_ax_noise = np.array([0.0, 5.0, 10.0, 20.0, 50.0, 100.0])
y_ax_SR = np.array([1.0, np.max(np.imag(perf_ring_eigvals_noise_5) / (-gamma/2)) / max_decay, np.max(np.imag(perf_ring_eigvals_noise_10) / (-gamma/2)) / max_decay, np.max(np.imag(perf_ring_eigvals_noise_20) / (-gamma/2)) / max_decay, np.max(np.imag(perf_ring_eigvals_noise_50) / (-gamma/2)) / max_decay, np.max(np.imag(perf_ring_eigvals_noise_100) / (-gamma/2)) / max_decay]) 

y_ax_10_noise = seg_10 / seg_10[0]
y_ax_20_noise = seg_20 / seg_20[0]
y_ax_50_noise = seg_50 / seg_50[0]
y_ax_80_noise = seg_80 / seg_80[0]

plt.scatter(x_ax_noise, y_ax_SR, c="purple")
plt.plot(x_ax_noise, y_ax_SR, c="purple")

plt.scatter(x_ax_noise, y_ax_80_noise, c="blue")
plt.plot(x_ax_noise, y_ax_80_noise, c="blue")

plt.scatter(x_ax_noise, y_ax_50_noise, c="green")
plt.plot(x_ax_noise, y_ax_50_noise, c="green")

plt.scatter(x_ax_noise, y_ax_20_noise, c="orange")
plt.plot(x_ax_noise, y_ax_20_noise, c="orange")

plt.scatter(x_ax_noise, y_ax_10_noise, c="red")
plt.plot(x_ax_noise, y_ax_10_noise, c="red")

plt.show()


