import sys
sys.path.append("..")
import data_tools as dt
from matplotlib import pyplot as plt
import numpy as np

mat = dt.matFromBed(sys.argv[1])

n = len(mat)

tots = np.zeros(n-1)
counts = np.zeros_like(tots)

for i in range(n):
	for j in range(i):
		s = i - j
		if mat[i,j] != 0:
			tots[s-1] += mat[i,j]
			counts[s-1] += 1

avgs = np.zeros_like(tots)

for i, (tot, count), in enumerate(zip(tots, counts)):
	if count != 0:
		avgs[i] = tot/count

plt.plot(range(n-1), avgs)
plt.xlabel("Separation (number of bins)")
plt.ylabel("Average contact frequency")
plt.show()
