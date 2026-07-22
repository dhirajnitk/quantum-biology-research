#!/bin/bash

MAX_DECAY_RATE_DIR="$1"
NOISE_VALUES="$2"
DECAY_RATE_FILE_PREFIX="$3"
POS_DIPOLE_FILE_PASSED_IN=$4
QUANTUM_YIELD_DIR="$5"

#
# Aggregate the results for all the noise runs by
# creating a numpy matrix with each row representing a specific noise
# value and the entire matrix normalized to the value of the case noise = 0.
#
NOISE_FILES_DECAY_RATE=""
for noise in ${NOISE_VALUES};
do
    if [[ $noise != "0" ]]; then
        NOISE_FILES_DECAY_RATE="${NOISE_FILES_DECAY_RATE} ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${noise}.log"
    fi
done

NOISE_FILES_QUANTUM_YIELD=""
for noise in ${NOISE_VALUES};
do
    NOISE_FILES_QUANTUM_YIELD="${NOISE_FILES_QUANTUM_YIELD} ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${noise}.log"
done

echo "$NOISE_FILES_DECAY_RATE"
echo
echo "$NOISE_FILES_QUANTUM_YIELD"

OUTPUT_FILENAME=${POS_DIPOLE_FILE_PASSED_IN}
echo "Position dipole file passed in: ${POS_DIPOLE_FILE_PASSED_IN}"
python create_plot_data.py \
    --decay-rate-files ${NOISE_FILES_DECAY_RATE} \
    --quantum-yield-files ${NOISE_FILES_QUANTUM_YIELD} \
    --zero-noise-file "${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}0.log" \
    --output-filename "${MAX_DECAY_RATE_DIR}/${OUTPUT_FILENAME}" \
    --noise_values ${NOISE_VALUES}
