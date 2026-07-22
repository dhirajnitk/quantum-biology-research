import argparse
import sys
import numpy as np

# Converts a list of quantum yield files (from input) to the quantum yield plot data.
def quantum_yield_files_to_plot_data(quantum_yield_files):
    plot_data_quantum_yield = [[]] * (NUM_NOISE_VALUES + 1)

    for quantum_yield_file in quantum_yield_files:
        with open(quantum_yield_file, "r") as fp:
            quantum_yield_contents = fp.readlines()
        quantum_yields_for_noise_value = []
        noise_value = None
        for each_line in quantum_yield_contents:
            quantum_yield = float(each_line.split(" ")[5].rstrip())
            # Verify spiral number
            noise_value = int(each_line.split(" ")[2])
            print(noise_value)
            quantum_yields_for_noise_value.append(quantum_yield)
        plot_data_quantum_yield[w_to_ind_qy[noise_value]] = quantum_yields_for_noise_value
        print()
        
    return plot_data_quantum_yield

# Converts a list of max decay rate files (from input) to the max decay rate plot data.
def max_decay_rate_files_to_plot_data(max_decay_rate_files):
    plot_data = [[]] * NUM_NOISE_VALUES
    normalized_plot_data = [[]] * NUM_NOISE_VALUES

    for max_decay_file in max_decay_rate_files:
        with open(max_decay_file, "r") as fp:
            max_decay_contents = fp.readlines()
        decay_rates_for_noise_value = []
        normalized_decay_rates_for_noise_value = []
        noise_value = None
        for each_line in max_decay_contents:
            decay_rate = float(each_line.split(" ")[5].rstrip())
            normalized_decay_rate = float(each_line.split(" ")[5].rstrip()) / no_noise_decay_rate
            # Verify spiral number
            noise_value = int(each_line.split(" ")[2])
            decay_rates_for_noise_value.append(decay_rate)
            normalized_decay_rates_for_noise_value.append(normalized_decay_rate)
        plot_data[w_to_ind[noise_value]] = decay_rates_for_noise_value
        normalized_plot_data[w_to_ind[noise_value]] = normalized_decay_rates_for_noise_value

    return plot_data, normalized_plot_data


# Parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument("--decay-rate-files", nargs="*", required=True, help="List of decay rate files with the results of multiple noise runs", )
parser.add_argument("--quantum-yield-files", nargs="*", required=True, help="List of quantum yield files with the results of multiple noise runs")
parser.add_argument("--zero-noise-file", required=True, help="File with contents of a single run with noise 0", )
parser.add_argument("--output-filename", required=True, help="Name of output file to write to")
parser.add_argument("--noise_values", nargs='+', required=True, help="The noise values")
args = parser.parse_args()

print(args.decay_rate_files)
print(args.zero_noise_file)
print(args.output_filename)
print(args.noise_values)

# Read the zero noise file
with open(args.zero_noise_file, "r") as fp:
    no_noise_contents = fp.readlines()
if len(no_noise_contents) != 1:
    print("ERROR: Zero noise file has [len(no_noise_contents)] lines but should only have 1 line")
    sys.exit(1)

# The file is structured like this: "5 spirals. 0 noise. Max_decay 31.80613136291504"
# Parse and get the spiral number and the decay rate for noise 0
no_noise_decay_rate = float(no_noise_contents[0].split(" ")[5].rstrip())
spiral_num = no_noise_contents[0].split(" ")[0]

# Now populate two dictionaries of noise values associated with indices (the one for QY including 0)
w_to_ind = {int(noise): index for index, noise in enumerate(args.noise_values[1:])}
w_to_ind_qy = {int(noise): index for index, noise in enumerate(args.noise_values)}

NUM_NOISE_VALUES = len(w_to_ind)

# If noise is zero i.e no noise
if NUM_NOISE_VALUES == 0:
    decay_rate_for_zero_noise = []
    decay_rate_for_zero_noise.append(no_noise_decay_rate)
    plot_data.append(decay_rate_for_zero_noise)
    normalized_plot_data.append(decay_rate_for_zero_noise)
else:
    # Generate max decay rate data
    plot_data, normalized_plot_data = max_decay_rate_files_to_plot_data(args.decay_rate_files)
    
    # Generate quantum yield data
    plot_data_quantum_yield = quantum_yield_files_to_plot_data(args.quantum_yield_files)

for elem in plot_data_quantum_yield:
    print(elem)

# Convert the list to a numpy array and save to a file.
np.save(f"{args.output_filename}_{spiral_num}_Spir.npy", np.array(plot_data))
np.save(f"{args.output_filename}-normalized_{spiral_num}_Spir.npy", np.array(normalized_plot_data))
np.save(f"{args.output_filename}-quantum-yield_{spiral_num}_Spir.npy", np.array(plot_data_quantum_yield, dtype=object))
