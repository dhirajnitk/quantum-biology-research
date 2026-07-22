import multiprocessing as mp
import numpy as np
import time
import sys

#np.set_printoptions(threshold=sys.maxsize)

def worker(arr_flattened):
	for i in range(100):
		for j in range(i+1):
			arr[i][j] = 100*i + j
			if i != j:
				arr[j][i] = arr[i][j]
	return arr

def mp_populate(arr):
	num_processes = mp.cpu_count()
	chunk_size = int((100*100) / num_processes)

	chunks = 

a = np.zeros((100, 100))
a = worker(a)

print(a)
print(mp.cpu_count())

