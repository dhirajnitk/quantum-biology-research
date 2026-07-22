#!/bin/bash

MAX_DECAY_RATE_DIR="$1"
NOISE="$2"
DECAY_RATE_FILE_PREFIX="$3"
QUANTUM_YIELD_DIR="$4"
#
# Combine the log files for the individual runs for a specific noise
#
echo "Noise = ${NOISE}"

rm -f ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}.log
rm -f ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}.log

# Decay Rate
noise_log_files=""
for i in $(seq 1 $NUM_REALIZATIONS); do
    noise_log_files="${noise_log_files} ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-{$i}.log"
done

echo "----"
echo "Max decay rate specific"
echo $noise_log_files
noise_log_files=$(ls ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-*.log)
echo "Max decay rate all"
echo $noise_log_files

for file in $noise_log_files; do 
    cat $file >> ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}.log
    echo "" >> ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}.log
done

# Quantum Yield
quantum_yield_log_files=""
for i in $(seq 1 $NUM_REALIZATIONS); do
    quantum_yield_log_files="${noise_log_files} ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-{$i}.log"
done

echo "QY specific"
echo $quantum_yield_log_files
quantum_yield_log_files=$(ls ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-*.log)
echo "QY all"
echo $quantum_yield_log_files
echo "----"

for file in $quantum_yield_log_files; do 
    cat $file >> ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}.log
    echo "" >> ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}.log
done
