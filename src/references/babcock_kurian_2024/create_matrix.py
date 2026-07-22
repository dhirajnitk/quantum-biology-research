"""
CREATED BY: HAMZA PATWA

This script creates the Hamiltonian matrix for a table containing positions and dipoles based on Celardo et. al. 2019.
Run the following command for usage instructions:
$ python create_matrix.py --help

"""

import argparse
import os.path
import sys
import time
from math import pi, sin, cos

import numpy as np
from numpy.linalg import norm

np.set_printoptions(threshold=sys.maxsize)

# Uncomment this line if you don't want numpy to print in scientific notation
# np.set_printoptions(suppress=True)


# CONSTANTS (in units described by footnotes 11 and 12 in Celardo et. al. 2019, New Journal of Physics).
h = 6.62607015e-27                # cm^2 g s^-1. Exact value of Planck constant from NIST
c = 2.99792458e10                 # cm s^-1.     Exact value of speed of light in vacuum from NIST
e0 = (1/280.0) * (1e7)            # cm^-1.       Trp excitation energy of 280 nm (converted to cm^-1 by taking the inverse and multiplying by 10^7)
mu = 6                            # Debeye (D).  Transition dipole strength of Trp.
k0 = 2*pi*e0 * (1e-8)             # A^-1.        Trp exitation wavenumber derived from e0
w0 = c*k0*(1e8)                   # s^-1.        Angular frequency of exitation of Trp derived from k0.
musq = ((1e-12)/(h*c))*(mu**2)    # A^3 cm^-1.   Square of the transition dipole strength of Trp (converted to diff units)
gamma = (4/3)*musq*(k0**3)        # cm^-1        Decay rate (0.0027304235183809407)
gamma_nr = 0.0183                 # cm^-1        Non-radiative decay rate
kB = 0.6950348                    # cm^-1 K^-1   Exact value of Boltzmann constant from NIST
T = 298                           # K            Room temperature
ns = 104

prg_start = time.perf_counter()

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", required=True, help="Filename of positions and dipoles file")
parser.add_argument("-d", required=True, help="Directory to store the output .npy file")
args = parser.parse_args()

# Extracting positions and dipole orientations from table
try:
    with open(args.f) as f:
        coordinates = f.readlines()

except FileNotFoundError:
    print(f"A file named {args.f} does not exist.")
    sys.exit(0)

num_dipoles = len(coordinates)
positions = np.array([[float(num) for num in line.split()[:3]] for line in coordinates])
dipoles = np.array([[float(num) for num in line.split()[3:]] for line in coordinates])


# Convenience functions to make it easier to construct matrix elements
def r(n, m):
    return norm(positions[n] - positions[m])


def r_hat(n, m):
    r_vec = positions[n] - positions[m]
    return r_vec / norm(r_vec)


def mu_hat(n):
    return dipoles[n]


def H0(n, m):
    if n == m:
        return (h * w0) / (2 * pi)
    return 0


def Delta(n, m):
    if n == m:
        return 0

    dimension_full_factor = (3 / 4) * gamma
    k0rnm = k0 * r(n, m)
    term1 = (
            ((-cos(k0rnm) / k0rnm) + (sin(k0rnm) / (k0rnm ** 2)) + (cos(k0rnm) / (k0rnm ** 3))) * mu_hat(n) @ mu_hat(m))
    term2 = ((-cos(k0rnm) / k0rnm) + (3 * sin(k0rnm) / (k0rnm ** 2)) + (3 * cos(k0rnm) / (k0rnm ** 3))) * (
            (mu_hat(n) @ r_hat(n, m)) * (mu_hat(m) @ r_hat(n, m)))
    return dimension_full_factor * (term1 - term2)


def G(n, m):
    if n == m:
        return gamma

    dimension_full_factor = (3 / 2) * gamma
    k0rnm = k0 * r(n, m)
    term1 = (((sin(k0rnm) / k0rnm) + (cos(k0rnm) / (k0rnm ** 2)) - (sin(k0rnm) / (k0rnm ** 3))) * mu_hat(n) @ mu_hat(m))
    term2 = ((sin(k0rnm) / k0rnm) + (3 * cos(k0rnm) / (k0rnm ** 2)) - (3 * sin(k0rnm) / (k0rnm ** 3))) * (
            (mu_hat(n) @ r_hat(n, m)) * (mu_hat(m) @ r_hat(n, m)))
    return dimension_full_factor * (term1 - term2)


# Populate the non-Hermitian Hamiltonian Heff using the above functions
pop_start = time.perf_counter()

Heff = np.zeros((num_dipoles, num_dipoles), dtype=np.complex64)

for n in range(num_dipoles):
    for m in range(n + 1):
        matrix_element = H0(n, m) + Delta(n, m) - (1j / 2) * G(n, m)
        Heff[n][m] = matrix_element
        if m != n:
            Heff[m][n] = matrix_element

fname_to_save = args.f.split("/")[-1].split(".")[0]

filepath_to_save = os.path.join(args.d, fname_to_save)
np.save(filepath_to_save, Heff)

pop_end = time.perf_counter()
print(f"Population time: {round(pop_end - pop_start, 3)}")
