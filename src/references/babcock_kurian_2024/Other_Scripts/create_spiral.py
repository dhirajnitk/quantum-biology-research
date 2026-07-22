'''
CREATED BY: HAMZA PATWA

This script creates all the .pdb files for a microtubule of a specified number of spirals.
Uses the "1jff.pdb" Protein DataBase (PDB) file.
Based on instructions in N. S. Babcock et. al. (2023) [Preprint] (arXiv:2302.01469) and G. L. Celardo et. al. (2019) (DOI: 10.1088/1367-2630/aaf839).

Run the following command for usage instructions:
$ python create_spiral.py --help 

'''

import numpy as np
from scipy.spatial.transform import Rotation as R
import math
import time
import sys
import argparse

'''
The read_pdb() method reads a PDB file.
-----------------------------
INPUT:  
	filename: Name of the PDB file as a string. 

OUTPUT: 
	(jffArr, jffAtoms): Tuple of length 2. The element at index 0 is an array that contains all the lines of the PDB file. The element at index 1 contains only those lines that start with the 'ATOM' keyword.

'''

def read_pdb(filename):
	
    with open(filename, "r") as jff:
        jffArr = jff.readlines()
	
    jffAtoms = [line for line in jffArr if line[0:4] == "ATOM"]
    return (jffArr, jffAtoms)


'''
The transform_dimer() method applies specified translations and rotations to one dimer (one PDB file).
-----------------------------
INPUT:  
	jffAtoms: Element at index 1 of the output of read_pdb().
	rotation_deg: The rotation angle in degrees. 
	translation: The translation amount in Angstroms.
	
OUTPUT: 
	jffAtomPositions: Array whose length is the number of atoms in the PDB file. Every element of the array contains another array of length 3, containing new x, y, and z positions of an atom in the PDB file.

'''

def transform_dimer(jffAtoms, rotation_deg, translation):
	
    # Each row of the below array contains the original (x, y, z) positions of an atom in "1jff.pdb."
    jffAtomPositions = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in jffAtoms])

    # Rotation by 55.38 degrees
    initial_rot = R.from_euler('x', 55.38, degrees=True)
    jffAtomPositions = initial_rot.apply(jffAtomPositions)

    # Rotation of -11.7 degrees about the beta tubulin Trp 346 CD2 atom (atom 5900 in the PDB file "1jff.pdb").
    initial_CD2_pos = jffAtomPositions[5898]
    jffAtomPositions = jffAtomPositions - initial_CD2_pos
    rot_around_CD2 = R.from_euler('x', -11.7, degrees=True)
    jffAtomPositions = rot_around_CD2.apply(jffAtomPositions)
    jffAtomPositions = jffAtomPositions + initial_CD2_pos

    # Translation by 112 Angstroms in the y direction and 3 Angstroms in the z direction.
    jffAtomPositions = jffAtomPositions + np.array([0.000, 112.000, 0.000])	
    jffAtomPositions = jffAtomPositions + np.array([0.000, 0.000, 3.000])

    # Extra translation to try and match with Nathan Babcock's version of Table A1 from Celardo et. al. 2019.
    # Uncomment the below line if you desire your table to match exactly with Nathan's table.
    # TODO: Maybe this line WOULD affect the physics, since it would lead to a microtubule with a different radius?
    #jffAtomPositions = jffAtomPositions + np.array([0.000, 0.959, -2.803])

    # Perform rotation and translation specified by input arguments
    r = R.from_euler('x', rotation_deg, degrees=True)
    jffAtomPositions = r.apply(jffAtomPositions)		
    jffAtomPositions = jffAtomPositions + np.array([translation, 0.000, 0.000])

    return jffAtomPositions
	

'''
The create_pdb() method creates a new PDB file from the new positions.
-----------------------------
INPUT:  
	originalArr: Element at index 0 of the output of read_pdb() (this is for all the other lines in the PDB file besides the ones starting with 'ATOM').
	originalArrAtoms: Element at index 1 of the output of read_pdb() (this is for the other characters in the lines starting with 'ATOM').
	newPositions: Output of transform_dimer(). These are substituted for the old positions.
	filename: Desired filename of the new PDB file that will be created.

OUTPUT: 
	NONE. A new PDB file is created, but this function doesn't return anything.

'''

def create_pdb(originalArr, originalArrAtoms, newPositions, filename):
	
    # TODO: This is very dangerous. PDB files should be parsed by the column number of the entry; values are NOT guarunteed to be separated by whitespace. 
    # TODO: The split() command is therefore risky.
    originalArrAtomsSplit = [line.split() for line in originalArrAtoms]

    # Array that will contain each line of the new PDB file, including the atom positions
    newArr = []

    # Array that will contain the new atom positions, with each element as a string.
    newArrAtoms = []

    # Counting index. Do not remove.
    counter = 0

    # Construct strings that will be the ones starting with 'ATOM'.
    for line, pos in zip(originalArrAtomsSplit, newPositions):
        rn = [min(8 - (len(str(pos[i]).split(".")[0])+1), 3) for i in range(3)]
        newArrAtoms.append(f"{line[0]}{line[1]:>7}  {line[2]:<4}{line[3]} {line[4]}{line[5]:>4}{pos[0]:>12.{rn[0]}f}{pos[1]:>8.{rn[1]}f}{pos[2]:>8.{rn[2]}f}  {line[9]}{line[10]:>6}           {line[11]}\n")

    # Construct all file lines of new PDB file.
    for line in originalArr:
        if line[0:4] != "ATOM":
            newArr.append(line)
        
        else:
            newArr.append(newArrAtoms[counter])
            counter = counter + 1

    # Create and write lines to the new PDB file.
    with open(filename, 'x') as newf:
        newf.writelines(newArr)
        
	
def main():
	
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, required=True, help="Number of spirals")
    args = parser.parse_args()
    numspirals = args.n

    """
    # Check if command line argument is present.
    if len(sys.argv) == 1:
        print("Please enter the number of spirals as a command line argument.")
        sys.exit(0)

    # start_time to time program. Can be commented if desired.
    start_time = time.time()
    """
    # Read the PDB file.
    jffArr, jffAtoms = read_pdb("PDB_files/1jff.pdb")
    """
    # Check to see if command line argument is an integer.
    try:
        numspirals = int(sys.argv[1])

    except ValueError:
        print(f"The argument {sys.argv[1]} is not an integer. Please enter an integer number of spirals.")
    """
    # Generate microtubule
    for j in range(numspirals):
        start_pos = 80*j
        for i in range(13):
            newPositions = transform_dimer(jffAtoms, -27.69*i, start_pos + 9*i)
            create_pdb(jffArr, jffAtoms, newPositions, f"1jff_num{(13*j)+i}.pdb")

    # Program timing.
    end_time = time.time()
    #print(f"Time for {numspirals} spirals: {round(end_time - start_time, 3)} seconds.")
	

if __name__ == "__main__":
	main()


