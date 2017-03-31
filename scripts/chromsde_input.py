import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

in_path = sys.argv[1]
mat_path = sys.argv[2]
id_path = sys.argv[3]

cluster = dt.clusterFromBed(in_path, None, None)
contactMat = dt.matFromBed(in_path, cluster)
n = len(contactMat)
maxNumDigits = int(np.ceil(np.log10(np.amax(contactMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(mat_path, contactMat, formatstring, delimiter="\t")

name = cluster.chrom.name
name = name[3:len(name)]	#remove "chr"
with open(id_path, "w") as out:
	for i, point in enumerate(cluster.getPoints()):
		out.write("\t".join((name, str(cluster.chrom.minPos + point.num*cluster.chrom.res), str(cluster.chrom.minPos + (point.num+1)*cluster.chrom.res), str(i+1))) + "\n")
out.close()
