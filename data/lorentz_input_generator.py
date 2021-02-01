import random
import numpy as np

n = 5000
d = 10

for i in range(n):
	data = np.random.randn(d-1)
	h = np.linalg.norm(data)
	h = np.sqrt(1 + h*h)
	print('%.10f\t' % h, end="")
	for j in range(d-1):
		print('%.10f' % data[j], end="")
		if j < d-2:
			print("\t", end="")
	if i < n-1:
		print("")