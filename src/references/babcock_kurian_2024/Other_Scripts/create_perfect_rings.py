'''
CREATED BY: HAMZA PATWA

This script creates a file with the positions and dipole orientations of a perfect ring structure with a set number of spirals.
The filename of the created file will be perfect_rings_table_pos_dip_{N}_spirals.txt, where N is the number of rings

Run the following command for usage instructions:
$ python create_perfect_rings.py --help 
'''

import numpy as np
from scipy.spatial.transform import Rotation as R
import math
import os
import sys
from matplotlib import pyplot as plt
import argparse

# Returns the midpoint between two points
def midpoint(point1, point2):
	return (point1 + point2)/2


# Returns the angle between two vectors
def angle_degrees(vec1, vec2):
	return np.arccos((vec1@vec2)/(math.sqrt(vec1@vec1) * math.sqrt(vec2@vec2))) * (180/math.pi)

# Returns a normalized version of the input vector
def normalize_vec(vec):
	return vec / math.sqrt(vec@vec)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, required=True, help='Number of rings')
    args = parser.parse_args()
    """
	if len(sys.argv) == 1:
		print("Please enter number of rings as a command line argument.")
		sys.exit(0)
	"""
	numrings = args.n
	dipoles_per_ring = 104
	angle = (2*math.pi) / dipoles_per_ring
	radius = 224 # units of angstroms
	
	positions = np.zeros((numrings * dipoles_per_ring, 3))
	
	# Paralell dipoles (PD) model
	dipole = np.array([0.0, 0.0, 1.0])
	dipoles = np.zeros((numrings * dipoles_per_ring, 3))

	for i in range(numrings):
		for j in range(dipoles_per_ring):
			if i == 0:
				dipoles[j,:] = normalize_vec(np.cross(np.array([0.0, 0.0, 1.0]), np.array([radius * math.cos(angle*j), radius * math.sin(angle*j), 0])))
			#dipoles[dipoles_per_ring*i + j,:] = dipoles[j,:]
			dipoles[dipoles_per_ring*i + j,:] = dipole
			positions[dipoles_per_ring*i + j,:] = np.array([radius * math.cos(angle*j), radius * math.sin(angle*j), i*80.0])
	'''
	# Tangent dipoles (TD) model
	dipole = np.array([0.0, 0.0, 1.0])
	dipoles = np.zeros((numrings * dipoles_per_ring, 3))

	for i in range(numrings):
		for j in range(dipoles_per_ring):
			if i == 0:
				dipoles[j,:] = normalize_vec(np.cross(np.array([0.0, 0.0, 1.0]), np.array([radius * math.cos(angle*j), radius * math.sin(angle*j), 0])))
			#dipoles[dipoles_per_ring*i + j,:] = dipoles[j,:]
			dipoles[dipoles_per_ring*i + j,:] = dipole
			positions[dipoles_per_ring*i + j,:] = np.array([radius * math.cos(angle*j), radius * math.sin(angle*j), i*80.0])
	'''
	
	# Create strings and write to file
	with open(f"perfect_rings_table_pos_dip_{args.n}_spirals.txt", 'x') as f:
		for dip, pos in zip(dipoles, positions):
			atoms_strings = f"{pos[0]:.3f} {pos[1]:.3f} {pos[2]:.3f} {dip[0]:.3f} {dip[1]:.3f} {dip[2]:.3f}\n"
			
			f.writelines(atoms_strings)
	
	'''
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	ax.quiver(positions[:,0], positions[:,1], positions[:,2], 30*dipoles[:,0], 30*dipoles[:,1], 30*dipoles[:,2], color='r')
	ax.scatter(positions[:,0], positions[:,1], positions[:,2], s=50)
	plt.show()
	'''
	
if __name__ == "__main__":
	main()


