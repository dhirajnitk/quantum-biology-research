#!/bin/bash

MAX_DECAY_RATE_DIR="$1"
LOCAL_NPY_MATRICES_DIRECTORY="$2"
POS_DIPOLE_FILENAME="$3"
NPY_MATRIX_FILE="$4"
QUANTUM_YIELD_DIR="$5"

echo "Create the Matrix if needed"

mkdir -p ${LOCAL_NPY_MATRICES_DIRECTORY}
mkdir -p ${MAX_DECAY_RATE_DIR}
mkdir -p ${QUANTUM_YIELD_DIR}

if [[ ! -f ${NPY_MATRIX_FILE}.npy ]];then
    echo "Matrix not found. Creating matrix..."
    python create_matrix.py -f ${POS_DIPOLE_FILENAME} -d ${LOCAL_NPY_MATRICES_DIRECTORY}
else
    echo "Matrix found - ${NPY_MATRIX_FILE}.npy"
    echo "Skipping the matrix creation step..."
fi
