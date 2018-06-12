import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

inpath = sys.argv[1]
outpath = sys.argv[2]

structure = dt.structureFromBed(inpath)
contactMat = dt.matFromBed(inpath, structure)
n = len(contactMat)
fullMat = np.zeros((n, n+2))

#locus IDs
for i, pointNum in enumerate(structure.getPointNums()):
	fullMat[i,0] = structure.chrom.minPos + structure.chrom.res * pointNum
	fullMat[i,1] = structure.chrom.minPos + structure.chrom.res * (pointNum + 1)

fullMat[:,2:n+2] = contactMat

maxNumDigits = int(np.ceil(np.log10(np.amax(fullMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(outpath, fullMat, formatstring, delimiter="\t")
