#!/bin/bash

#
# This is specific to the arc setup.
#
module load anaconda3
python3 --version ; pip3 --version
pip3 install numpy matplotlib torch
