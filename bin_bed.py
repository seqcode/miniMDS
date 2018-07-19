import data_tools as dt
import sys
import numpy as np

path = sys.argv[1]
res = int(sys.argv[2])
outpath = sys.argv[3]

chrom = dt.chromFromBed(path)
chrom.res = res
chrom.minPos = int(np.floor(float(chrom.minPos)/res)) * res	#round
chrom.maxPos = int(np.ceil(float(chrom.maxPos)/res)) * res

struct = dt.structureFromBed(path, chrom)
mat = dt.matFromBed(path, struct)

points = struct.getPoints()

with open(outpath, "w") as out:
	for i in range(len(mat)):
		abs_index1 = points[i].absolute_index
		for j in range(i):
			if mat[i,j] != 0:
				abs_index2 = points[j].absolute_index
				out.write("\t".join((chrom.name, str(chrom.getGenCoord(abs_index1)), str(chrom.getGenCoord(abs_index1) + res), chrom.name, str(chrom.getGenCoord(abs_index2)), str(chrom.getGenCoord(abs_index2) + res), str(mat[i,j]))))
				out.write("\n")
	out.close()
