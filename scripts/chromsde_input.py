import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

in_path = sys.argv[1]
mat_path = sys.argv[2]
id_path = sys.argv[3]

structure = dt.structureFromBed(in_path, None, None)
contactMat = dt.matFromBed(in_path, structure)
n = len(contactMat)
maxNumDigits = int(np.ceil(np.log10(np.amax(contactMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(mat_path, contactMat, formatstring, delimiter="\t")

name = structure.chrom.name
name = name[3:len(name)]	#remove "chr"
with open(id_path, "w") as out:
	for i, point in enumerate(structure.getPoints()):
		out.write("\t".join((name, str(structure.chrom.minPos + point.num*structure.chrom.res), str(structure.chrom.minPos + (point.num+1)*structure.chrom.res), str(i+1))) + "\n")
out.close()
