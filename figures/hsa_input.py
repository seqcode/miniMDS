import numpy as np
import sys
sys.path.append("../../minimds")
from tools import Tracker

inpath = sys.argv[1]
outpath = sys.argv[2]

def matFromBed(path):
	"""Create array from BED"""
	minPos = sys.float_info.max
	maxPos = 0
	count = 0
	res = None
	print "Scanning {}".format(path)
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			pos1 = int(line[1])
			pos2 = int(line[4])
			if res is None:
				res = int(line[2]) - pos1
			if pos1 < minPos:
				minPos = pos1
			if pos2 > maxPos:
				maxPos = pos2
			count += 1
	infile.close()
	
	n = (maxPos - minPos)/res + 1
	mat = np.zeros((n, n+2))

	for i in range(n):
		mat[i,0] = minPos + res*i
		mat[i,1] = minPos + res*i+1

	tracker = Tracker("Reading {}".format(path), count)
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			pos1 = int(line[1])
			pos2 = int(line[4])
			binNum1 = (pos1 - minPos)/res
			binNum2 = (pos2 - minPos)/res
			if binNum1 is not None and binNum2 is not None:
				value = float(line[6])
				mat[binNum1, binNum2 + 2] += value
				mat[binNum2, binNum1 + 2] += value
			tracker.increment()
	infile.close()
	return mat

contactMat = matFromBed(inpath)
maxNumDigits = int(np.ceil(np.log10(np.amax(contactMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(outpath, contactMat, formatstring, delimiter="\t")
