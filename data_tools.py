import sys
import numpy as np
from tools import Tracker
import linear_algebra as la
import array_tools as at

class ChromParameters(object):
	def __init__(self, minPos, maxPos, res, name, size):
		self.minPos = minPos	#minimum genomic coordinate
		self.maxPos = maxPos	#maximum genomic coordinate
		self.res = res		#resolution (bp)
		self.name = name	#e.g. "chr22"
		self.size = size	#number of lines in file

class Cluster(object):
	"""Intrachromosomal cluster of points or subclusters in 3-D space"""
	def __init__(self, points, clusters, chrom, offset):
		self.points = points
		self.clusters = clusters	#subclusters
		for cluster in self.clusters:	#auto-fill
			for point in cluster.points:
				self.points.append(point)	
		self.chrom = chrom
		self.offset = offset	#indexing offset (for subclusters only)

	def getCoords(self):
		return [point.pos for point in self.getPoints()]

	def getPointNums(self):
		return np.array([point.num for point in self.getPoints()])

	def getPoints(self):
		return self.points[np.where(self.points !=0)[0]]

	def getIndex(self, genCoord):
		"""Converts genomic coordinate into index"""
		pointNum = self.chrom.getPointNum(genCoord)
		if pointNum is None:
			return None
		else:
			pointNum = pointNum - self.offset
			if pointNum >= 0 and pointNum < len(self.points):
				point = self.points[pointNum]
				if point == 0:
					return None
				else:
					return point.index
			else:
				return None
	
	def setClusters(self, clusters):
		self.clusters = clusters
		self.points = np.zeros(max([max(cluster.getPointNums()) for cluster in clusters]) + 1, dtype=np.object)	#reset
		for cluster in self.clusters:
			for point in cluster.points:
				if point != 0:
					self.points[point.num] = point

	def createSubcluster(self, points, offset):
		"""Creates subcluster containing pointsToAdd"""
		subcluster = Cluster(points, [], self.chrom, offset)
		self.clusters.append(subcluster)

	def transform(self, r, t, reflect=False):
		"""Rotates by r; translates by t; reflect if True"""
		if r is None:	#default: no rotation
			r = np.mat(np.identity(3))
		if t is None:	#default: no translation
			t = np.mat(np.zeros(3)).T
		a = np.mat(self.getCoords())
		n = len(a)
		a_transformed = np.array(((r*a.T) + np.tile(t, (1, n))).T)
		for i, pointNum in enumerate(self.getPointNums()):
			self.points[pointNum - self.offset].pos = a_transformed[i]
		if reflect:
			self.reflection(0)

	def reflection(self, axis):
		"""Reflects cluster across given axis"""
		reflectionVector = np.ones(3)
		reflectionVector[axis] = -1
		for point in self.getPoints():
			point.pos = np.multiply(reflectionVector, point.pos)

	def write(self, outpath):
		with open(outpath, "w") as out:
			out.write(self.chrom.name + "\n")
			out.write(str(self.chrom.res) + "\n")
			out.write(str(self.chrom.minPos) + "\n")
			num = self.offset
			for point in self.points:
				if point == 0:
					out.write("\t".join((str(num), "nan", "nan", "nan")) + "\n")
				else:
					out.write("\t".join((str(num), str(point.pos[0]), str(point.pos[1]), str(point.pos[2]))) + "\n")
				num += 1
		out.close()
	def indexPoints(self):
		i = 0
		for point in self.points:
			if point != 0:
				point.index = i
				i += 1

class Point(object):
	"""Point in 3-D space"""
	def __init__(self, pos, num, chrom, index):
		self.pos = pos	#3D coordinates
		self.num = num	#locus (not necessarily sequential)
		self.chrom = chrom	#chromosome parameters
		self.index = index	#sequential

def clusterFromBed(path, chrom, tads):
	"""Initializes cluster from intrachromosomal BED file."""
	if chrom is None:
		chrom = intraChromFromBed(path, None)

	cluster = Cluster([], [], chrom, 0)
	
	#get TAD for every locus
	if tads is None:
		tadNums = np.zeros(cluster.chrom.getLength())
	else:
		tadNums = []
		tadNum = 1
		for tad in tads:
			for i in range(tad[0], tad[1]):
				tadNums.append(tadNum)
			tadNum += 1
	maxIndex = len(tadNums) - 1

	points_to_add = np.zeros(cluster.chrom.getLength(), dtype=np.bool)	#true if locus should be added
	tracker = Tracker("Identifying loci", cluster.chrom.size)

	#find which loci should be added
	with open(path) as listFile:
		for line in listFile:
			line = line.strip().split()
			pos1 = int(line[1])
			pos2 = int(line[4])
			pointNum1 = cluster.chrom.getPointNum(pos1)
			pointNum2 = cluster.chrom.getPointNum(pos2)
			if pointNum1 is not None and pointNum2 is not None:
				tadNum1 = tadNums[min(pointNum1, maxIndex)]
				tadNum2 = tadNums[min(pointNum2, maxIndex)]
				if pointNum1 != pointNum2 and tadNum1 == tadNum2:		#must be in same TAD
					if points_to_add[pointNum1] == False:
						points_to_add[pointNum1] = True
					if points_to_add[pointNum2] == False:
						points_to_add[pointNum2] = True
			tracker.increment()
	listFile.close()

	#create points
	points = np.zeros(cluster.chrom.getLength(), dtype=np.object)
	pointNums = np.where(points_to_add == True)[0]
	for pointNum in pointNums:
		points[pointNum] = Point((0,0,0), pointNum, cluster.chrom, None)
	cluster.points = points
	cluster.indexPoints()
	
	return cluster

def intraChromFromBed(path, res):
	"""Initialize ChromParams from intrachromosomal file in BED format"""
	count = 0
	minPos = sys.float_info.max
	maxPos = 0
	print "Scanning {}".format(path)
	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			if count == 0:
				name = line[0]
			pos1 = int(line[1])
			pos2 = int(line[4])
			if pos1 < minPos:
				minPos = pos1
			elif pos1 > maxPos:
				maxPos = pos1
			if pos2 < minPos:
				minPos = pos2
			elif pos2 > maxPos:
				maxPos = pos2
			if res is None:
				res = (int(line[2]) - pos1)	
			count += 1
	infile.close()
	minPos = minPos/res * res	#round
	maxPos = maxPos/res * res
	return ChromParameters(minPos, maxPos, res, name, count)

def basicParamsFromBed(path):
	size = 0
	res = None
	print "Scanning {}".format(path)
	with open(path) as infile:
		for line in infile:
			size += 1
			if res is None:
				line = line.strip().split()
				res = (int(line[2]) - int(line[1]))
	infile.close()
	return size, res

def matFromBed(path, cluster):	
	"""Converts BED file to matrix. Only includes loci in cluster."""
	cluster.indexPoints()
	pointNums = cluster.getPointNums()

	numpoints = len(pointNums)
	maxPointNum = max(pointNums)
	minPointNum = min(pointNums)
	mat = np.zeros((numpoints, numpoints))	

	assert maxPointNum - cluster.offset < len(cluster.points)

	with open(path) as infile:
		for line in infile:
			linearray = line.strip().split()	#line as array of strings
			loc1 = int(linearray[1])
			loc2 = int(linearray[4])
			index1 = cluster.getIndex(loc1)
			index2 = cluster.getIndex(loc2)
			if index1 is not None and index2 is not None:
				if index1 > index2:
					row = index1
					col = index2
				else:
					row = index2
					col = index1
				mat[row, col] += float(linearray[6])
	infile.close()


	at.makeSymmetric(mat)
	rowsums = np.array([sum(row) for row in mat])
	empty = np.where(rowsums == 0)[0]
	assert len(np.where(rowsums == 0)[0]) == 0

	at.makeSymmetric(mat)

	return mat

def highToLow(highCluster, resRatio):
	"""Reduces resolution of cluster"""
	lowChrom = highCluster.chrom.reduceRes(resRatio)

	low_n = int(np.ceil(len(highCluster.points)/float(resRatio)))

	lowCluster = Cluster(np.zeros(low_n, dtype=np.object), [], lowChrom, highCluster.offset/resRatio)

	allPointsToMerge = []
	for i in range(len(lowCluster.points)):
		allPointsToMerge.append([])
	
	for highPoint in highCluster.getPoints():
		pointsToMerge = []
		highNum = highPoint.num
		lowNum = highNum/resRatio
		allPointsToMerge[lowNum - lowCluster.offset].append(highPoint)

	index = lowCluster.offset
	for i, pointsToMerge in enumerate(allPointsToMerge):
		if len(pointsToMerge) > 0:
			lowCluster.points[i] = mergePoints(pointsToMerge, i + lowCluster.offset, lowChrom, index)
			index += 1

	return lowCluster

def mergePoints(pointsToMerge, newPointNum, chrom, index):
	"""Creates new point with average position of pointsToMerge"""
	coords = np.array([point.pos for point in pointsToMerge])
	meanCoord = np.mean(coords, axis=0)
	return Point(meanCoord, newPointNum, chrom, index)

def clusterFromFile(path):
	hasMore = True
	with open(path) as infile:
		name = infile.readline().strip()
		res = int(infile.readline().strip())
		minPos = int(infile.readline().strip())
		chrom = ChromParameters(minPos, None, res, name, None)
		cluster = Cluster([], [], chrom, 0)
		index = 0
		while hasMore:
			line = infile.readline().strip().split()
			if len(line) == 0:
				hasMore = False
			else:
				num = int(line[0])
				if line[1] == "nan":
					point = 0
				else:
					x = float(line[1])
					y = float(line[2])
					z = float(line[3])
					point = Point((x,y,z), num, chrom, index)
					index += 1
				cluster.points.append(point)
	infile.close()
	cluster.points = np.array(cluster.points)
	cluster.chrom.maxPos = cluster.chrom.minPos + cluster.chrom.res*num	#max pos is last point num
	return cluster
