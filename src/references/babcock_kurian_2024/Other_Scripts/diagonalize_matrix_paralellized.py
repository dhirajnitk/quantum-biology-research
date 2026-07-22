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
import multiprocessing
import multiprocessing.managers

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

# Setup the custom Manager that allows different subprocesses to access shared NumPy array
class MyManager(multiprocessing.managers.BaseManager):
    pass

MyManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)

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
		return (h*w0) / (2*pi)
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

# Function to parallelize
def populate_row_Heff(n):
	
	# Utilize the fact that Heff is symmetric (we don't need to evaluate the same matrix element twice)
	for i in range(n+1):
		matrix_element = H0(n,i) + Delta(n,i) - (1j/2)*G(n,i)
		Heff[n,i] = matrix_element
		
		if i != n:
			Heff[i,n] = matrix_element


#def populate_row_Heff(n):
#	for i in range(n+1):
#		elem = num_dipoles*n + i
#		Heff[n,i] = elem
#
#		if n != i:
#			Heff[i,n] = elem

# MAIN METHOD WOULD START HERE

m = MyManager()
m.start()

# Initialize array with Manager
Heff = m.np_zeros((num_dipoles, num_dipoles), dtype=np.complex64)

# Initialize processes
pool = multiprocessing.Pool()

# Multiprocesssing!
mp_start = time.perf_counter()
#processes = [pool.apply_async(populate_row_Heff, args=(i,)) for i in range(num_dipoles)]
#results = [p.get() for p in processes]
results = pool.map(populate_row_Heff, range(num_dipoles))
mp_end = time.perf_counter()

# Now the array should be populated! The rest of the program continues as normal.

# Diagonalize matrix
diag_start = time.time()
eigenvalues, eigenvectors = eig(Heff)
diag_end = time.time()

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
print(f"Population time for {int(num_dipoles/104)} spirals: {round(mp_end - mp_start, 3)}.")

#print(index_of_SRS(eigenvalues))
#print()
plt.xlim(-150, 150)
plt.ylim(-0.5, 20)

plt.xlabel("Energy")
plt.ylabel("Γ / γ")
	
plt.scatter(energy_values, Gamma_values)
plt.show()
	


