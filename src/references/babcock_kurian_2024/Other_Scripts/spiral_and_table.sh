#!/bin/bash

# This script runs create_spiral.py and create_table.py concurrently, cleaning up afterwards.
# This script takes in 2 command line arguments:
# The first is the number of spirals and the second is the output filename (the argumnents to create_spiral.py and create_table.py, respectively)

numspirals=$1
pos_dip_fname=$2

python create_spiral.py -n $numspirals
python create_table.py -f $pos_dip_fname
rm 1jff_num*

