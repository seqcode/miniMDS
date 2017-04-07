import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

inpath = sys.argv[1]
outpath = sys.argv[2]

cluster = dt.clusterFromBed(inpath, None, None)
contactMat = dt.matFromBed(inpath, cluster)
n = len(contactMat)
fullMat = np.zeros((n, n+2))

#locus IDs
for i, pointNum in enumerate(cluster.getPointNums()):
	fullMat[i,0] = cluster.chrom.minPos + cluster.chrom.res * pointNum
	fullMat[i,1] = cluster.chrom.minPos + cluster.chrom.res * (pointNum + 1)

fullMat[:,2:n] = contactMat

maxNumDigits = int(np.ceil(np.log10(np.amax(contactMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(outpath, contactMat, formatstring, delimiter="\t")
