import random
import numpy as np

n = 20
d = 10

for i in range(n):
	data = np.random.randn(d)
	data = data / np.linalg.norm(data) * np.random.rand()  # norm must be within 1
	for j in range(d):
		print('%.10f' % data[j], end="")
		if j < d-1:
			print("\t", end="")
	if i < n-1:
		print("")