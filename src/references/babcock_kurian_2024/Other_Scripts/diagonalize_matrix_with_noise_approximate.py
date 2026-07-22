'''
CREATED BY: HAMZA PATWA

This script diagonalizes the effective Hamiltonian Heff from Celardo et. al. 2019. It takes in a command line argument that is the filename
of the file that contains dipole positions and orientations. It should be created with create_table.py BEFORE running this script.

USAGE: $ python diagonalize_matrix.py filename

filename is the name of the file that contains dipole positions and orientations.
'''

from numpy.linalg import eig, eigvals, norm, inv
import numpy as np
from math import pi, sin, cos
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import sys

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
W = 0                      # cm^-1        (energy)
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
		
	
if len(sys.argv) == 1:
	print("Please enter the name of the file with dipole positions and orientations as a command-line argument.")
	sys.exit(0)

prg_start = time.time()

# Extracting positions and dipole orientations from table
try:
	with open(sys.argv[1]) as f:
		coordinates = f.readlines()

except FileNotFoundError:
	print(f"A file named {sys.argv[1]} does not exist.")
	sys.exit(0)

num_dipoles = len(coordinates)
num_spirals = int(num_dipoles/ns)
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
		return np.random.uniform(H0elem - W, H0elem + W)
		
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

def j_proj_q(j, q):
	if j not in range(1, 105) or q not in range(1, 105):
		print("ERROR")
		sys.exit(0)
	
	return np.exp((1j*2*np.pi*j*q)/104)

def s_proj_r(s, r):
	if s not in range(1, num_spirals+1) or r not in range(1, num_spirals+1):
		print("ERROR")
		sys.exit(0)

	return np.sin((np.pi*s*r) / (num_spirals+1))


# Populate the non-Hermitian Hamiltonian Heff using the above
Heff = np.zeros((num_dipoles, num_dipoles), dtype=np.complex64)

# Store diagonal energies to maybe use again
noiseW = np.zeros(num_dipoles, dtype=np.complex64)

# Utilize the fact that Heff is symmetric (we don't need to evaluate the same matrix element twice)
for n in range(num_dipoles):
	for m in range(n+1):
		site_e = H0(n,m)
		matrix_element = site_e + Delta(n,m) - (1j/2)*G(n,m)
		Heff[n][m] = matrix_element
		if m != n:
			Heff[m][n] = matrix_element
	
	noiseW[n] = site_e


# Diagonalize matrix
diag_start = time.time()

U = np.zeros((num_dipoles, num_dipoles), dtype=np.complex64)

fact_eig = (1/np.sqrt(104)) * (np.sqrt(2)/np.sqrt(num_spirals + 1))
ring_proj_q1 = np.array([j_proj_q(j, 1) for j in range(1, 105)])

for r in range(1, num_spirals+1):
    for q in range(1, 105):

        ring_proj_q = ring_proj_q1**q
        chain_proj_r = np.array([s_proj_r(s, r) for s in range(1, num_spirals+1)])

        eigst = fact_eig * np.kron(chain_proj_r, ring_proj_q)
        U[:,104*(r-1) + (q-1)] = eigst


#print(U[:10,0] / ((1/np.sqrt(104)) * (np.sqrt(2)/np.sqrt(num_spirals + 1))))

approx_diag_Heff = inv(U)@Heff@U

diag_end = time.time()

eigenvalues = np.diag(approx_diag_Heff)
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
plt.xlim(-150, 150)
plt.ylim(-0.5, 20)

plt.xlabel("Energy")
plt.ylabel("Γ / γ")
	
plt.scatter(energy_values, Gamma_values)
plt.show()
	


