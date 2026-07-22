"""
CREATED BY: HAMZA PATWA

This script diagonalizes the effective Hamiltonian Heff from Celardo et. al. 2019 with thermal noise.
It takes in the Hamiltonian created by create_matrix.py.
It creates a .log file that contains a single value: the maximum decay rate.

USAGE: $ python diagonalize_matrix.py filename.npy

filename is the name of the Hamiltonian .npy file.
"""

import argparse
import sys
import time
import os
from math import pi

import numpy as np
import torch

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"GPU: {torch.cuda.is_available()}")
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
gamma = (4/3)*musq*(k0**3)        # cm^-1        Decay rate of Trp (0.0027304235183809407)
gamma_nr = (gamma/0.124)-gamma    # cm^-1        Non-radiative decay rate of Trp (0.019289120984691162). Using the experimental QY of Trp in solution of 0.124.
kB = 0.6950348                    # cm^-1 K^-1   Exact value of Boltzmann constant from NIST
T = 298                           # K            Room temperature

# I tried to convert everything to SI units just as a sanity check (I should get the same results) but I didn't get
# the same results (I'm probably doing something wrong here).

"""
h = 6.626e-13                # cm^2 g s^-1
c = 2.998e17                 # cm s^-1
e0 = 7.094e2                 # nm^2 g s^-2  (energy)
k0 = 2.24e-2                 # nm^-1        (wavenumber)
gamma = 5.423e-5             # nm^2 g s^-2  (energy)
w0 = c*k0                    # s^-1         (angular frequency)
ns = 104
"""

parser = argparse.ArgumentParser()
parser.add_argument("-W", type=int, required=True, help="Integer value of noise")
parser.add_argument("-H", required=True, help="Name of matrix (Hamiltonian) file to read from")
parser.add_argument("-dro", required=True, help="Decay Rate Output (DRO): name of output file to write max decay rate to")
parser.add_argument("-qyo", required=True, help="Quantum Yield Output (QYO): name of output file to write quantum yield to")
parser.add_argument("--drfp", required=True, help="Name of the environment variable DECAY_RATE_FILE_PREFIX")
parser.add_argument("--qyfp", required=True, help="Name of the environment variable QUANTUM_YIELD_FILE_PREFIX")
parser.add_argument("--trp_per_spiral", required=True, help="Number of Trp molecules per spiral")
parser.add_argument("--dir_to_save_eigvals", help="Name of directory to save eigenvalues, leave blank if you don't want to save")

args = parser.parse_args()
W = args.W
hamiltonian_npy = args.H

prg_start = time.time()
print(f"noise: {W}")


# Returns the location of the superradiant (SR) state in the spectrum of energy eigenvalues as an index (integer).
# 1 means that the SR state is the ground state, and 104 (for one spiral) would mean the SR state is the most excited
#  state.
# The expression with np.where() that's being returned is just a complicated-looking way to pick out the index.


def print_time(t):
    hrs = np.floor(t / 3600.0)
    mins = np.floor(((t / 3600.0) - hrs) * 60)
    secs = np.round((((t / 3600.0) - hrs) * 60 - mins) * 60, 2)

    if hrs == 0:
        if mins == 0:
            return f"{secs} sec."
        else:
            return f"{mins} min, {secs} sec."
    else:
        return f"{hrs} hr, {mins} min, {secs} sec."


def index_of_SRS(eigenvalues):
    Energ_eigenvalues = np.real(eigenvalues)
    Gamma_eigenvalues = np.imag(eigenvalues) / (-gamma / 2)
    return np.where(np.sort(Energ_eigenvalues) == Energ_eigenvalues[
        np.where(Gamma_eigenvalues == np.max(Gamma_eigenvalues))[0][0]])[0][0]


def H0(n, m):
    if n == m:
        H0elem = (h * w0) / (2 * pi)
        # Add noise W
        return np.random.uniform(H0elem - (W / 2.0), H0elem + (W / 2.0))

    return 0


# Return a random value in the range [value - noise/2, value + noise/2]
def apply_noise(value, noise):
    return np.random.uniform(value - (noise/2.0), value + (noise/2.0))


# Calculates thermal quantum yield given complex eigenvalues
def thermal_QY(eigenvalues):
    Gamma_values = np.imag(eigenvalues) / (-1/2)
    energy_values = np.real(eigenvalues)
    
    Z = np.sum(np.array([np.exp(-Ej/(kB*T)) for Ej in energy_values]))
    thermal_avg_Gamma = np.sum(np.array([Gammaj * np.exp(-Ej/(kB*T)) for Gammaj, Ej in zip(Gamma_values, energy_values)])) / Z 
    thermal_avg_QY = thermal_avg_Gamma / (thermal_avg_Gamma + gamma_nr)
    return thermal_avg_QY


# Load matrix, and change diagonal elements
Heff = np.load(hamiltonian_npy)
numspirals = float(Heff.shape[0]) / int(args.trp_per_spiral)

for i in range(np.shape(Heff)[0]):
    Heff[i][i] = apply_noise(np.real(Heff[i][i]), W) + (np.imag(Heff[i][i]) * 1j)

print("start diagonalizing...")
diag_start = time.perf_counter()

# if torch.cuda.is_available():
#     Heff_tensor = torch.from_numpy(Heff).to(device=device)
#     eigenvalues = torch.linalg.eigvals(Heff_tensor).cuda()
# else:
#     Heff_tensor = torch.from_numpy(Heff).to(device=device)
#     eigenvalues = torch.linalg.eigvals(Heff_tensor)

# Heff_tensor = torch.from_numpy(Heff).to(device=device)
# eigenvalues = torch.linalg.eigvals(Heff_tensor)

eigenvalues = np.linalg.eigvals(Heff)

diag_end = time.perf_counter()
print("end diagonalizing.")

# Gamma_values = torch.imag(eigenvalues) / (-gamma / 2)
# energy_values = torch.real(eigenvalues)

Gamma_values = np.imag(eigenvalues) / (-gamma/2)
energy_values = np.real(eigenvalues)

prg_end = time.time()

print(f"Program runtime: {print_time(prg_end - prg_start)}.")
print(f"Diagonalization time for {numspirals} spirals: {print_time(diag_end - diag_start)}.")

with open(args.dro, "w") as f:
    f.write(f"{numspirals:.5f} spirals. {W} noise. Max_decay {np.max(Gamma_values)}")

# eigenvalues = eigenvalues.detach().cpu().numpy()
print(f"Sum of decay rates: {np.sum(Gamma_values)}")

with open(args.qyo, "w") as f:
    f.write(f"{numspirals:.5f} spirals. {W} noise. Quantum_yield {thermal_QY(eigenvalues)}")

if args.dir_to_save_eigvals is not None:
    if W == 0:
        np.save(f"{args.dir_to_save_eigvals}/{numspirals:.5f}_Spir_MT_eigvals", eigenvalues)

print("DONE")
