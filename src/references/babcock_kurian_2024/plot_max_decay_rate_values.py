import argparse
import sys
import distinctipy
import numpy as np
from matplotlib import pyplot as plt
import os

plt.rcParams["font.family"] = "Nimbus Roman"
plt.rcParams["font.size"] = 20

trp_to_structname = {
    24: "amyloid6MST",
    36: "amyloid6DSO",
    608: "actin6BNO_19_MT",
    104: "microtubules"
}

t = int(os.environ['NUM_SPIRALS_TO_DISCARD_IN_NOISE_ANALYSIS'])
trp_per_spiral = int(os.environ['TRP_PER_SPIRAL'])

def spiral_num_from_fname(fname):
    return round(float(fname.split("_")[-2]), 5)
    # spiral_num = 0
    # for part in fname.split("_"):
    #     if part.isdigit():
    #         spiral_num = int(part)
    #         break
    # if spiral_num == 0:
    #     print("ERROR: spiral number was not detected")
    #     sys.exit(-1)
    # return spiral_num


parser = argparse.ArgumentParser()
parser.add_argument("--decay_rate_files", nargs="+", help="Numpy array filenames that contain the maximum decay rate values for different noises.")
parser.add_argument("--quantum_yield_files", nargs="+", help="Numpy array filenames that contain the quantum yield values for different noises.")
parser.add_argument("--noise_values", nargs='+', type=float, required=True, help="The noise values")
args = parser.parse_args()

NUM_NOISE_VALUES = len(args.noise_values) - 1
NUM_SPIRALS = len(args.decay_rate_files)
print(f"NUM_NOISE_VALUES: {NUM_NOISE_VALUES}")
print(f"NUM_SPIRALS: {NUM_SPIRALS}")

# Extract spiral numbers and create dictionary that maps spiral number to index
spiral_numbers = [spiral_num_from_fname(fname) for fname in args.decay_rate_files]
sp_to_ind = {sp_num: i for i, sp_num in enumerate(sorted(spiral_numbers))}
ind_to_sp = {ind: sp_num for sp_num, ind in sp_to_ind.items()}

noise_values_excluding_0 = np.array(args.noise_values)[1:]
print(noise_values_excluding_0)

# Array of numpy files. Each numpy file contains all realizations for all noises for one spiral
all_max_decay_rate_files = [None] * NUM_SPIRALS
all_quantum_yield_files = [None] * NUM_SPIRALS

# 1D array of values. Contains the zero-noise quantum yield values for each spiral
quantum_yield_zero_noise_values = [None] * NUM_SPIRALS

# Loads numpy arrays into program
#print(args.decay_rate_files)
for fname in args.decay_rate_files:
    nam = fname.split("/")[-1]
    print(f"Name of numpy file --> {nam}")
    spiral_num = spiral_num_from_fname(fname)
    all_max_decay_rate_files[sp_to_ind[spiral_num]] = np.load(fname)

for fname in args.quantum_yield_files:
    nam = fname.split("/")[-1]
    print(f"Name of quantum yield file --> {nam}")
    spiral_num = spiral_num_from_fname(fname)
    np_loaded_arr = np.load(fname, allow_pickle=True)
    all_quantum_yield_files[sp_to_ind[spiral_num]] = np.array([l for l in np_loaded_arr[1:]])
    quantum_yield_zero_noise_values[sp_to_ind[spiral_num]] = np_loaded_arr[0][0]

max_decay_rate_averages = np.array([np.mean(all_realizations_for_one_spiral, axis=1) for all_realizations_for_one_spiral in all_max_decay_rate_files])
quantum_yield_averages = np.array(
    [np.insert(np.mean(all_realizations_for_one_spiral, axis=1), 0, zero_noise_quantum_yield) for (all_realizations_for_one_spiral, zero_noise_quantum_yield) in zip(all_quantum_yield_files, quantum_yield_zero_noise_values)]
)
errors = np.array([np.std(all_realizations_for_one_spiral, axis=1) for all_realizations_for_one_spiral in all_max_decay_rate_files])
errors_quantum_yield = np.array([np.insert(np.std(all_realizations_for_one_spiral, axis=1), 0, 0) for all_realizations_for_one_spiral in all_quantum_yield_files])

# Sets up and does figure
plt.figure(figsize=(12, 8))
plt.ylim(0, 1.1)
# Size of error bar
capsize_val = 7
colors = ["red", "orange", "yellow", "green", "blue", "purple", "brown", "gray", "black", "pink", "cyan", "magenta", "lightblue", "navy"]

for i in range(len(max_decay_rate_averages)):
    print(f"color for {i} is {colors[i]}")
    print(noise_values_excluding_0)
    print(max_decay_rate_averages[i])
    plt.scatter(noise_values_excluding_0, max_decay_rate_averages[i], c=colors[i], label=f"{ind_to_sp[i]} spirals")
    plt.errorbar(noise_values_excluding_0, max_decay_rate_averages[i], yerr=errors[i], capsize=capsize_val, c=colors[i])

plt.xlabel("Noise")
plt.ylabel("Max. Γ(W)/Γ(0)")
plt.legend()
plt.tight_layout()
plt.savefig(f"max_decay_rate_plot_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"max_decay_rate_plot_{trp_to_structname[trp_per_spiral]}.png")

plt.xscale('log')
plt.savefig(f"max_decay_rate_plot_log_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"max_decay_rate_plot_log_{trp_to_structname[trp_per_spiral]}.png")

# Quantum Yield plot (no zero noise tracing back)
plt.figure(figsize=(12, 8))

for i in range(NUM_NOISE_VALUES + 1):
    print(f"color for {i} is {colors[i]}")
    print(sorted(spiral_numbers))
    print(quantum_yield_averages[:,i])
    plt.scatter(np.array(sorted(spiral_numbers))[t:]*trp_per_spiral, quantum_yield_averages[:,i][t:], c=colors[i], label=f"W = {args.noise_values[i]}")
    plt.errorbar(np.array(sorted(spiral_numbers))[t:]*trp_per_spiral, quantum_yield_averages[:,i][t:], yerr=errors_quantum_yield[:,i][t:], capsize=capsize_val, c=colors[i])

plt.xlabel("# of Trp")
plt.ylabel("Quantum Yield")
plt.legend()
plt.tight_layout()
plt.savefig(f"quantum_yield_plot_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"quantum_yield_plot_{trp_to_structname[trp_per_spiral]}.png")

plt.xscale('log')
plt.savefig(f"quantum_yield_plot_log_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"quantum_yield_plot_log_{trp_to_structname[trp_per_spiral]}.png")

# Quantum yield plots for zero noise tracing back
plt.figure(figsize=(12, 8))

for i in range(NUM_NOISE_VALUES + 1):
    print(f"color for {i} is {colors[i]}")
    print(sorted(spiral_numbers))
    print(quantum_yield_averages[:,i])
    if args.noise_values[i] != 0:
        plt.scatter(np.array(sorted(spiral_numbers)[t:])*trp_per_spiral, quantum_yield_averages[:,i][t:], c=colors[i], label=f"W = {args.noise_values[i]}")
        plt.errorbar(np.array(sorted(spiral_numbers)[t:])*trp_per_spiral, quantum_yield_averages[:,i][t:], yerr=errors_quantum_yield[:,i][t:], capsize=capsize_val, c=colors[i])
    else:
        plt.scatter(np.array(sorted(spiral_numbers))*trp_per_spiral, quantum_yield_averages[:,i], c=colors[i], label=f"W = {args.noise_values[i]}")
        plt.errorbar(np.array(sorted(spiral_numbers))*trp_per_spiral, quantum_yield_averages[:,i], yerr=errors_quantum_yield[:,i], capsize=capsize_val, c=colors[i])

plt.xlabel("# of Trp")
plt.ylabel("Quantum Yield")
plt.legend()
plt.tight_layout()
plt.savefig(f"quantum_yield_plot_zero_noise_tracing_back_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"quantum_yield_plot_zero_noise_tracing_back_{trp_to_structname[trp_per_spiral]}.png")

plt.xscale('log')
plt.savefig(f"quantum_yield_plot_log_zero_noise_tracing_back_{trp_to_structname[trp_per_spiral]}.pdf")
plt.savefig(f"quantum_yield_plot_log_zero_noise_tracing_back_{trp_to_structname[trp_per_spiral]}.png")