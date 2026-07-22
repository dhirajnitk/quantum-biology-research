import numpy as np


with open("PAPER_table_positions_dipoles_wrong.txt", 'r') as f:
	filelines = f.readlines()
	
filesplit = [[float(num) for num in line.split()] for line in filelines]
#index 4,6,8
modified_file = [f"{fs[0]} {fs[1]} {fs[2]} {fs[3]} {fs[5]} {fs[7]}\n" for fs in filesplit]

with open("Positions_Dipoles_Tables/PAPER_table_positions_dipoles.txt", 'x') as f:
	f.writelines(modified_file)

