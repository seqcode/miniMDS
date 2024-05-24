import sys
sys.path.append("..")
import data_tools as dt
import numpy as np

path = sys.argv[1]
res = int(sys.argv[2])
outpath = sys.argv[3]

chrom1, chrom2 = dt.chromFromBed(path, True)

chrom1.res = res
chrom1.minPos = dt.round_down(chrom1.minPos, res)
chrom1.maxPos = dt.round_up(chrom1.maxPos, res)
struct1 = dt.structureFromBed(path, chrom=chrom1, chrom_order=1)
points1 = struct1.getPoints()

chrom2.res = res
chrom2.minPos = dt.round_down(chrom2.minPos, res)
chrom2.maxPos = dt.round_up(chrom2.maxPos, res)
struct2 = dt.structureFromBed(path, chrom=chrom2, chrom_order=2)
points2 = struct2.getPoints()

mat = dt.matFromBed(path, structure1=struct1, structure2=struct2)

with open(outpath, "w") as out:
	for i in range(mat.shape[0]):
		gen_coord1 = chrom1.getGenCoord(points1[i].absolute_index)
		for j in range(mat.shape[1]):
			if mat[i,j] != 0:
				gen_coord2 = chrom2.getGenCoord(points2[j].absolute_index)
				out.write("\t".join((chrom1.name, str(gen_coord1), str(gen_coord1 + res), chrom2.name, str(gen_coord2), str(gen_coord2 + res), str(mat[i,j]))))
				out.write("\n")