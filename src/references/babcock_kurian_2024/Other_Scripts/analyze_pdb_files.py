# Location 3228 has a TER

def main():
	
	# Creating arrays with all the atoms
	with open("PDB_files/1jff.pdb", "r") as jff:
		jffArr = jff.readlines()
	jffAtoms = [line for line in jffArr if line[0:4] == "ATOM"]
	
	with open("PDB_files/1tub.pdb", "r") as tub:
		tubArr = tub.readlines()
	tubAtoms = [line for line in tubArr if line[0:4] == "ATOM"]

	# Total number of atoms in both of them (matches with ChimeraX)
	print(len(jffAtoms) + len(tubAtoms))
	
	# Total number of each amino acid in each
	amino_acids = [*set([atom.split()[3] for atom in jffAtoms])]
	amino_acids.sort()
	
	jff_amino_acid_counts = [0]*20
		
	for atom in jffAtoms:
		amino_acid_index = amino_acids.index(atom.split()[3])
		jff_amino_acid_counts[amino_acid_index] = jff_amino_acid_counts[amino_acid_index] + 1
	
	tub_amino_acid_counts = [0]*20
		
	for atom in tubAtoms:
		amino_acid_index = amino_acids.index(atom.split()[3])
		tub_amino_acid_counts[amino_acid_index] = tub_amino_acid_counts[amino_acid_index] + 1

	print("1JFF")	
	print(amino_acids)	
	print(jff_amino_acid_counts)
	
	print()
	
	print("1TUB")
	print(amino_acids)	
	print(tub_amino_acid_counts)
	
	print()
	
	print("DIFF")
	print([b-a for a, b in zip(jff_amino_acid_counts, tub_amino_acid_counts)])
	
if __name__ == "__main__":
	main()
