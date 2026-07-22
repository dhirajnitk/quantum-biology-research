import numpy as np
import math

def normalize(vec):
	return vec / math.sqrt(vec@vec)

def distance(vec):
	return math.sqrt(sum([val**2 for val in vec]))

def angle(vec1, vec2):
	return (180/math.pi) * np.arccos(vec1@vec2 / (math.sqrt(vec1@vec1)*math.sqrt(vec2@vec2)))

def main():
	with open("hamza_table_A1_from_preprint_changed_sign_of_rotation.txt", 'r') as f:
		my_pos = f.readlines()
	
	with open("Positions_Dipoles_Tables/PAPER_table_positions_dipoles.txt", 'r') as f:
		paper_pos = f.readlines()
	
	with open("nathan_table.txt", 'r') as f:
		nathan_pos = f.readlines()
	
	my_pos = np.array([[float(num) for num in pos.split()] for pos in my_pos])
	paper_pos = np.array([[float(num) for num in pos.split()] for pos in paper_pos])
	nathan_pos = np.array([[float(num) for num in pos.split()] for pos in nathan_pos])

	#print(my_pos[0:8])
	
	for i in range(8):
		d = np.round(my_pos[i] - paper_pos[np.where(paper_pos[:,0] == my_pos[i][0])[0][0]], 2)
		m = my_pos[i] 
		p = paper_pos[np.where(paper_pos[:,0] == my_pos[i][0])[0][0]]
		dn = np.round(my_pos[i] - nathan_pos[i], 5)

		
		#print(f"{m[0]:<10}{m[1]:<10}{m[2]:<10}{m[3]:<10}{m[4]:<10}{m[5]:<10}")
		#print(f"{p[0]:<10}{p[1]:<10}{p[2]:<10}{p[3]:<10}{p[4]:<10}{p[5]:<10}")
		print(f"{dn[0]:<10}{dn[1]:<10}{dn[2]:<10}{dn[3]:<10}{dn[4]:<10}{dn[5]:<10}")
		
		#print(f"{d[0]:<10}{d[1]:<10}{d[2]}")
		#print(f"{round(angle(m,p), 2)} deg")
		print()
		#print(f"{d[0]:<10}{d[1]:<10}{d[2]}")
		#print(distance(my_pos[i]))
		#print(distance(paper_pos[np.where(paper_pos[:,0] == my_pos[i][0])[0][0]]))
		
	
	print()
	
	my_pos = my_pos[np.argsort(my_pos[:,0])]
	paper_pos = paper_pos[np.argsort(paper_pos[:,0])]
		
	#print(my_pos[0:5,:])
	#print(paper_pos)
	#print()
	#print(np.round(my_pos - paper_pos, 2))	
	#print()
	#print((180/math.pi) * my_pos[0]@paper_pos[0] / (math.sqrt(my_pos[0]@my_pos[0])*math.sqrt(paper_pos[0]@paper_pos[0])))
	#print([round(distance(pos1) - distance(pos2), 2) for pos1, pos2 in zip(my_pos, paper_pos)])

if __name__ == "__main__":
	main()
