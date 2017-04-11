import sys
sys.path.append("..")
import data_tools as dt
from tools import Tracker
import simple_tad as tad
import heatmap as hm
import numpy as np
import misc

def matFromDixon(path, chrom):
	"""Creates contact matrix from Dixon tsv data"""
	numBins = misc.getLength(chrom)
	mat = np.zeros((numBins, numBins))
	tracker = Tracker("Reading " + path, chrom.size)
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			pos1 = int(line[0])
			pos2 = int(line[1])
			if pos1 != pos2:
				if pos1 >= chrom.minPos and pos1 <= chrom.maxPos and pos2 >= chrom.minPos and pos2 <= chrom.maxPos:
					bin1 = misc.getPointNum(chrom, pos1)	
					bin2 = misc.getPointNum(chrom, pos2)
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

def plotDixon(mat):
	tads = [[0,8], [8,38], [38,52], [52,78], [78,97], [97,115], [115,127]]
	outpath = "Fig1A"
	hm.heatMapFromMat(mat, 100, tads, outpath)

def plotMovingAverage(mat):
	smoothingFactor = 5
	outpath = "Fig1B"
	tads = tad.getDomains(mat, smoothingFactor, 0)
	hm.heatMapFromMat(mat, 100, tads, outpath) 

minPos = 49000000	#from Dixon
maxPos = 54066692	#from Dixon
res = 40000	#from Dixon
name = "chr22"
size = 30949158
path = "mESC_chr6.tsv"

chrom = dt.ChromParameters(minPos, maxPos, res, name, size)

mat = matFromDixon(path, chrom)
plotDixon(mat)
plotMovingAverage(mat)
