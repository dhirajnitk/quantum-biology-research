'''
CREATED BY: HAMZA PATWA

This script diagonalizes the effective Hamiltonian Heff from Celardo et. al. 2019. It takes in a command line argument that is the filename
of the file that contains dipole positions and orientations. It should be created with create_table.py BEFORE running this script.

Run the following command for usage instructions:
$ python diagonalize_pre_created_matrix_with_noise.py --help
'''

from numpy.linalg import eig, eigvals, norm
import numpy as np
from math import pi, sin, cos
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import sys
import argparse
import torch

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"GPU: {torch.cuda.is_available()}")

np.set_printoptions(threshold=sys.maxsize)

# Uncomment this line if you don't want numpy to print in scientific notation
#np.set_printoptions(suppress=True)

# CONSTANTS (in units described by footnotes 11 and 12 in Celardo et. al. 2019, New Journal of Physics).

h = 6.626e-27                # cm^2 g s^-1
c = 2.998e10                 # cm s^-1
e0 = 35716.65                # cm^-1        (energy)
k0 = 2.24e-3                 # A^-1         (inverse length)
musq = 181224                # A^3 cm^-1 
gamma = 2.73e-3              # cm^-1        (energy)
w0 = c*k0*(1e8)              # s^-1         (angular frequency)
ns = 104
#W = 5                      # cm^-1        (energy)
#num_dipoles = ns*int(sys.argv[1])

# I tried to convert everything to SI units just as a sanity check (I should get the same results) but I didn't get the same results (I'm probably doing something wrong here).

'''
h = 6.626e-13                # cm^2 g s^-1
c = 2.998e17                 # cm s^-1
e0 = 7.094e2                 # nm^2 g s^-2  (energy)
k0 = 2.24e-2                 # nm^-1        (wavenumber)
gamma = 5.423e-5             # nm^2 g s^-2  (energy)
w0 = c*k0                    # s^-1         (angular frequency)
ns = 104
'''

parser = argparse.ArgumentParser()
parser.add_argument("-W", type=int, required=True, help="Value of thermal noise")
args = parser.parse_args()
W = args.W

print(f"noise: {W}")

# Returns the location of the superradiant (SR) state in the spectrum of energy eigenvalues as an index (integer).
# 1 means that the SR state is the ground state, and 104 (for one spiral) would mean the SR state is the most excited state.
# The expression with np.where() that's being returned is just a complicated-looking way to pick out the index.

def index_of_SRS(eigenvalues):
    Energ_eigenvalues = np.real(eigenvalues)
    Gamma_eigenvalues = np.imag(eigenvalues) / (-gamma/2)
    return np.where(np.sort(Energ_eigenvalues) == Energ_eigenvalues[np.where(Gamma_eigenvalues == np.max(Gamma_eigenvalues))[0][0]])[0][0]

def H0(n,m):
	if n == m:
		H0elem = (h*w0) / (2*pi)
		# Add noise W
		return np.random.uniform(H0elem - (W/2.0), H0elem + (W/2.0))
		
	return 0

# Load matrix, and change diagonal elements
Heff = np.load("Eigenvalues/perf_ring_Hamiltonian_100_spirals_noise_0.npy")

for i in range(np.shape(Heff)[0]):
	Heff[i][i] = H0(i,i)

print("start diag...")
# Diagonalize matrix
diag_start = time.time()
#eigenvalues, eigenvectors = eig(Heff)
#eigenvalues = eigvals(Heff)

if torch.cuda.is_available():
	eigenvalues = torch.linalg.eigvals(torch.from_numpy(Heff).to(device=device)).cuda()
else:
	eigenvalues = torch.linalg.eigvals(torch.from_numpy(Heff).to(device=device))

diag_end = time.time()

print("diag finished.")
Gamma_values = np.imag(eigenvalues) / (-gamma/2)
energy_values = np.real(eigenvalues)

prg_end = time.time()

#SRS = eigenvectors[:,index_of_SRS(eigenvalues)]

#things_to_plot = np.concatenate((positions, np.abs(SRS[:,np.newaxis])**2), axis=1)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#scatterplot = ax.scatter(things_to_plot[:,0], things_to_plot[:,1], things_to_plot[:,2], c=things_to_plot[:,3], cmap='jet', vmin=1e-6, vmax=1e-2, s=200, alpha=1)

#colorbar_ax = fig.add_subplot(122)
#plt.colorbar(scatterplot, ax=ax)

#plt.show()
print(f"Program runtime: {round(prg_end - prg_start, 3)}.")
print(f"Diagonalization time for {int(num_dipoles/104)} spirals: {round(diag_end - diag_start, 3)}.")

#print(index_of_SRS(eigenvalues))
#print()

#np.save(f"perf_rings_100_spiral_noise_{W}_eigenvalues", eigenvalues)


plt.xlim(-150, 150)
plt.ylim(-0.5, 20)

plt.xlabel("Energy")
plt.ylabel("Γ / γ")
	
plt.scatter(energy_values, Gamma_values)
plt.show()
	

