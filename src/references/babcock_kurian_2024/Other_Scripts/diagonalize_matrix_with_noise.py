'''
CREATED BY: HAMZA PATWA

This script diagonalizes the effective Hamiltonian Heff from Celardo et. al. 2019. It takes in a command line argument that is the filename
of the file that contains dipole positions and orientations. It should be created with create_table.py BEFORE running this script.

USAGE: $ python diagonalize_matrix.py filename

filename is the name of the file that contains dipole positions and orientations.
'''

from numpy.linalg import eig, eigvals, norm
import numpy as np
from math import pi, sin, cos
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import sys
import argparse

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
W = 10                      # cm^-1        (energy)
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

# Returns the location of the superradiant (SR) state in the spectrum of energy eigenvalues as an index (integer).
# 1 means that the SR state is the ground state, and 104 (for one spiral) would mean the SR state is the most excited state.
# The expression with np.where() that's being returned is just a complicated-looking way to pick out the index.

def index_of_SRS(eigenvalues):
    Energ_eigenvalues = np.real(eigenvalues)
    Gamma_eigenvalues = np.imag(eigenvalues) / (-gamma/2)
    return np.where(np.sort(Energ_eigenvalues) == Energ_eigenvalues[np.where(Gamma_eigenvalues == np.max(Gamma_eigenvalues))[0][0]])[0][0]
		

prg_start = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("-W", type=int, required=True)
parser.add_argument("-f", required=True)
args = parser.parse_args()
W = args.W

# Extracting positions and dipole orientations from table
try:
	with open(args.f) as f:
		coordinates = f.readlines()

except FileNotFoundError:
	print(f"A file named {args.f} does not exist.")
	sys.exit(1)


num_dipoles = len(coordinates)
positions = np.array([[float(num) for num in line.split()[:3]] for line in coordinates])
dipoles = np.array([[float(num) for num in line.split()[3:]] for line in coordinates])
#dipoles = np.array([[0.0, 0.0, float(line.split()[5])] for line in coordinates])

# Convenience functions to make it easier to construct matrix elements
def r(n,m):
	return norm(positions[n] - positions[m])

def r_hat(n,m):
	r_vec = positions[n] - positions[m]
	return r_vec / norm(r_vec)

def mu_hat(n):
	return dipoles[n]

def H0(n,m):
	if n == m:
		H0elem = (h*w0) / (2*pi)
		# Add noise W
		return np.random.uniform(H0elem - (W/2.0), H0elem + (W/2.0))
		
	return 0
	
def Delta(n,m):		
	if n == m:
		return 0

	dimension_full_factor = (3/4)*gamma
	k0rnm = k0*r(n,m)
	term1 = ((-cos(k0rnm) / k0rnm) + (sin(k0rnm) / (k0rnm**2)) + (cos(k0rnm) / (k0rnm**3))) * mu_hat(n)@mu_hat(m)
	term2 = ((-cos(k0rnm) / k0rnm) + (3*sin(k0rnm) / (k0rnm**2)) + (3*cos(k0rnm) / (k0rnm**3))) * ((mu_hat(n)@r_hat(n,m)) * (mu_hat(m)@r_hat(n,m)))
	return dimension_full_factor * (term1 - term2)

	#return (mu_hat(n)@mu_hat(m) - 3*((mu_hat(n)@r_hat(n,m)) * (mu_hat(m)@r_hat(n,m)))) / (r(n,m)**3)

	
def G(n,m):
	if n == m:
		return gamma
	
	dimension_full_factor = (3/2)*gamma
	k0rnm = k0*r(n,m)
	term1 = ((sin(k0rnm) / k0rnm) + (cos(k0rnm) / (k0rnm**2)) - (sin(k0rnm) / (k0rnm**3))) * mu_hat(n)@mu_hat(m)
	term2 = ((sin(k0rnm) / k0rnm) + (3*cos(k0rnm) / (k0rnm**2)) - (3*sin(k0rnm) / (k0rnm**3))) * ((mu_hat(n)@r_hat(n,m)) * (mu_hat(m)@r_hat(n,m)))
	return dimension_full_factor * (term1 - term2)
	
	#return gamma * mu_hat(n)@mu_hat(m)


# Populate the non-Hermitian Hamiltonian Heff using the above
Heff = np.zeros((num_dipoles, num_dipoles), dtype=np.complex64)

print("Populating matrix...")
# Utilize the fact that Heff is symmetric (we don't need to evaluate the same matrix element twice)
for n in range(num_dipoles):
	for m in range(n+1):
		site_e = H0(n,m)
		matrix_element = site_e + Delta(n,m) - (1j/2)*G(n,m)
		Heff[n][m] = matrix_element
		if m != n:
			Heff[m][n] = matrix_element
	
print("Done populating. Now diagonalizing...")
# Diagonalize matrix
diag_start = time.time()
#eigenvalues, eigenvectors = eig(Heff)
eigenvalues = eigvals(Heff)
diag_end = time.time()
print("Done diagonalizing.\n")
Gamma_values = np.imag(eigenvalues) / (-gamma/2)
energy_values = np.real(eigenvalues)

prg_end = time.time()

print(f"Program runtime: {round(prg_end - prg_start, 3)}.")
print(f"Diagonalization time for {int(num_dipoles/104)} spirals: {round(diag_end - diag_start, 3)}.")

spiral_num = (args.f).split("_")[6]

np.save(f"perf_ring_eigenvalues_{spiral_num}_spirals_noise_{W}", eigenvalues)

if W == 0:
	np.save(f"perf_ring_Hamiltonian_{spiral_num}_spirals_noise_{W}", Heff)

