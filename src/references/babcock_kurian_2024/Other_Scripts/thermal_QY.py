from numpy import array, complex64
import numpy as np
from matplotlib import pyplot as plt
import os
import sys

h = 6.626e-27                # cm^2 g s^-1  (Planck constant)
c = 2.998e10                 # cm s^-1      (Speed of light in vacuum)
e0 = 35716.65                # cm^-1        (energy)
k0 = 2.24e-3                 # A^-1         (inverse length)
musq = 181224                # A^3 cm^-1
gamma = 2.73e-3              # cm^-1        (energy)
w0 = c*k0*(1e8)              # s^-1         (angular frequency)
kB = 0.695                   # cm^1 K^-1    (Boltzmann constant)
T = 298                      # K            (Room temperature)
ns = 104
gammanr = 0.0183             # cm^-1

AX_LAB_SIZE = 25
AX_TICK_SIZE = 20
AX_LEGEND_SIZE = 15

def index_of_SRS(eigenvalues):
	Energ_eigenvalues = np.real(eigenvalues)
	Gamma_eigenvalues = np.imag(eigenvalues) / (-gamma/2)
	return np.where(np.sort(Energ_eigenvalues) == Energ_eigenvalues[np.where(Gamma_eigenvalues == np.max(Gamma_eigenvalues))[0][0]])[0][0]


def thermal_QY(eigenvalues):
	Gamma_values = np.imag(eigenvalues) / (-gamma/2)
	energy_values = np.real(eigenvalues)
	
	Z = np.sum(np.array([np.exp(-Ej/(kB*T)) for Ej in energy_values]))
	return np.sum(np.array([Gammaj * np.exp(-Ej/(kB*T)) for Gammaj, Ej in zip(Gamma_values, energy_values)])) / Z
	

with open("Eigenvalues/perfect_rings_PD_eigenvalues_numerical_diag_100_spiral.txt", 'r') as f:
	eigenvals_PD_str = f.readlines()
	eigenvals_PD = eval(" ".join(eigenvals_PD_str).replace("\n", ""))

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


directory = "Positions_Dipoles_Tables/"
seg_10 = []
seg_20 = []
seg_50 = []
seg_80 = []

for filename in os.listdir(directory):
	f = os.path.join(directory, filename)
	if os.path.isfile(f):
		if f[-4:] == ".npy" and "eigenvalues" in f:
			print(f)


plt.figure(figsize=(12,8))
plt.scatter(np.real(eigenvals_PD), np.imag(eigenvals_PD) / (-gamma/2), c='b')
plt.title("PD")

print(f"PD (SR high): {thermal_QY(eigenvals_PD) / (thermal_QY(eigenvals_PD) + gammanr)}")

plt.figure(figsize=(12,8))

QYn = []
QYn.append(thermal_QY(perf_ring_eigvals_noise_5))
QYn.append(thermal_QY(perf_ring_eigvals_noise_10))
QYn.append(thermal_QY(perf_ring_eigvals_noise_20))
QYn.append(thermal_QY(perf_ring_eigvals_noise_50))
QYn.append(thermal_QY(perf_ring_eigvals_noise_100))

plt.scatter(range(5), QYn)

plt.show()







