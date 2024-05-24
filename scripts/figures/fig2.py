import sys
sys.path.append("..")
from data_tools import ChromParameters
from tools import Tracker
import heatmap as hm
import simple_tad as tad
import numpy as np

def matFromDixon(path, chrom):
	"""Creates contact matrix from Dixon tsv data"""
	numBins = chrom.getLength()
	mat = np.zeros((numBins, numBins))
	tracker = Tracker("Reading " + path, chrom.size)
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			pos1 = int(line[0])
			pos2 = int(line[1])
			if pos1 != pos2:
				if pos1 >= chrom.minPos and pos1 <= chrom.maxPos and pos2 >= chrom.minPos and pos2 <= chrom.maxPos:
					bin1 = chrom.getAbsoluteIndex(pos1)	
					bin2 = chrom.getAbsoluteIndex(pos2)
					if bin1 > bin2:
						row = bin1
						col = bin2
					else:
						row = bin1
						col = bin2
					mat[row, col] += 1
			tracker.increment()
	infile.close()
	return mat

def plotLevels(mat):
	smoothingFactors = [1, 2, 3, 8, 33]	#these smoothing factors were selected to demonstrate to best demonstrate TAD levels
	domainsToInclude = [range(1, 15), [2,3,4,5], [7], [1,6], [3]]	#selected domains from these smoothing factors to maximize prettiness	
	all_tads = []
	for i in range(len(smoothingFactors)):	
		smoothingFactor = smoothingFactors[i]
		indices = domainsToInclude[i]
		tads = tad.getDomains(mat, smoothingFactor, 0)
		for index in indices:
			all_tads.append(tads[index])
	hm.heatMapFromMat(mat, 100, all_tads, "Fig2")	#all levels combined

minPos = 49000000	#from Dixon
maxPos = 54066692	#from Dixon
res = 40000	#from Dixon
name = "chr22"
size = 30949158
path = "mESC_chr6.tsv"

chrom = ChromParameters(minPos, maxPos, res, name, size)

mat = matFromDixon(path, chrom)
plotLevels(mat)
