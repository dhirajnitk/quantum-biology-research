#!/bin/bash

MAX_DECAY_RATE_DIR="$1"
NPY_MATRIX_FILE="$2"
NOISE="$3"
CI_NODE_INDEX="$4"
DECAY_RATE_FILE_PREFIX="$5"
QUANTUM_YIELD_DIR="$6"
TRP_PER_SPIRAL="$7"
DIR_TO_SAVE_EIGVALS="$8"

#
# Check if the log file exists. If so then skip diagonalizing the matrix
#
if [ -f ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log ] && \
   [ -f ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log ]; then
    echo "The files ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log" 
    echo "and ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log"
    echo "already exists, so skip diagonalizing the matrix"
else
    echo "Diagonalizing Hamiltonian"
    echo "Noise = ${NOISE}"

    if [ -z "$DIR_TO_SAVE_EIGVALS" ]; then
        python diagonalize_pre_created_matrix_with_noise_returning_max_decay_rate.py \
        -W ${NOISE} \
        -H ${NPY_MATRIX_FILE}.npy \
        -dro ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log \
        -qyo ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log \
        --drfp ${DECAY_RATE_FILE_PREFIX} \
        --qyfp ${QUANTUM_YIELD_FILE_PREFIX} \
        --trp_per_spiral ${TRP_PER_SPIRAL}
    else
        python diagonalize_pre_created_matrix_with_noise_returning_max_decay_rate.py \
            -W ${NOISE} \
            -H ${NPY_MATRIX_FILE}.npy \
            -dro ${MAX_DECAY_RATE_DIR}/${DECAY_RATE_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log \
            -qyo ${QUANTUM_YIELD_DIR}/${QUANTUM_YIELD_FILE_PREFIX}${NOISE}-${CI_NODE_INDEX}.log \
            --drfp ${DECAY_RATE_FILE_PREFIX} \
            --qyfp ${QUANTUM_YIELD_FILE_PREFIX} \
            --trp_per_spiral ${TRP_PER_SPIRAL} \
            --dir_to_save_eigvals ${DIR_TO_SAVE_EIGVALS}
    fi
fi
