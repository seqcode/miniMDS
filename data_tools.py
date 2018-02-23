import sys
import numpy as np
from tools import Tracker
import linear_algebra as la
import array_tools as at

class ChromParameters(object):
	"""Basic information on chromosome, inferred from input file"""
	def __init__(self, minPos, maxPos, res, name, size):
		self.minPos = minPos	#minimum genomic coordinate
		self.maxPos = maxPos	#maximum genomic coordinate
		self.res = res		#resolution (bp)
		self.name = name	#e.g. "chr22"
		self.size = size	#number of lines in file

	def getLength(self):
		"""Number of possible loci"""
		return (self.maxPos - self.minPos)/self.res + 1

	def getPointNum(self, genCoord):
		"""Converts genomic coordinate into point number"""
		if genCoord < self.minPos or genCoord > self.maxPos:
			return None
		else:
			return int((genCoord - self.minPos)/self.res) 

	def reduceRes(self, resRatio):
		"""Creates low-res version of this chromosome"""
		lowRes = self.res * resRatio
		lowMinPos = (self.minPos/lowRes)*lowRes		#approximate at low resolution
		lowMaxPos = (self.maxPos/lowRes)*lowRes
		return ChromParameters(lowMinPos, lowMaxPos, lowRes, self.name, self.size)

class Structure(object):
	"""Intrachromosomal structure of points or substructures in 3-D space"""
	def __init__(self, points, structures, chrom, offset):
		self.points = points
		self.structures = structures	#substructures
		for structure in self.structures:	#auto-fill
			for point in structure.points:
				self.points.append(point)	
		self.chrom = chrom	#chromosome parameters
		self.offset = offset	#indexing offset (for substructures only)

	def getCoords(self):
		return [point.pos for point in self.getPoints()]

	def setCoords(self, coords):
		for coord, point_num in zip(coords, self.getPointNums()):
			self.points[point_num].pos = coord

	def getPointNums(self):
		return np.array([point.num for point in self.getPoints()])

	def getPoints(self):
		return self.points[np.where(self.points !=0)[0]]

	def getGenCoords(self):
		"""Non-null genomic coordinates of structure"""
		return [self.chrom.minPos + self.chrom.res * point_num for point_num in self.getPointNums()]

	def getIndex(self, genCoord):
		"""Converts genomic coordinate into index"""
		pointNum = self.chrom.getPointNum(genCoord)
		if pointNum is None:
			return None
		else:
			pointNum -= self.offset
			if pointNum >= 0 and pointNum < len(self.points):
				point = self.points[pointNum]
				if point == 0:
					return None
				else:
					return point.index
			else:
				return None
	
	def setstructures(self, structures):
		self.structures = structures
		self.points = np.zeros(max([max(structure.getPointNums()) for structure in structures]) + 1, dtype=np.object)	#reset
		for structure in self.structures:
			for point in structure.points:
				if point != 0:
					self.points[point.num] = point

	def createSubstructure(self, points, offset):
		"""Creates substructure containing points"""
		substructure = Structure(points, [], self.chrom, offset)
		self.structures.append(substructure)

	def transform(self, r, t):
		"""Rotates by r; translates by t"""
		if r is None:	#default: no rotation
			r = np.mat(np.identity(3))
		if t is None:	#default: no translation
			t = np.mat(np.zeros(3)).T
		a = np.mat(self.getCoords())
		n = len(a)
		a_transformed = np.array(((r*a.T) + np.tile(t, (1, n))).T)
		for i, pointNum in enumerate(self.getPointNums()):
			self.points[pointNum - self.offset].pos = a_transformed[i]

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

class Point(object):
	"""Point in 3-D space"""
	def __init__(self, pos, num, chrom, index):
		self.pos = pos	#3D coordinates
		self.num = num	#locus (not necessarily sequential)
		self.chrom = chrom	#chromosome parameters
		self.index = index	#sequential

def structureFromBed(path, chrom, tads):
	"""Initializes structure from intrachromosomal BED file."""
	if chrom is None:
		chrom = chromFromBed(path)

	structure = Structure([], [], chrom, 0)
	
	#get TAD for every locus
	if tads is None:
		tadNums = np.zeros(structure.chrom.getLength())
	else:
		tadNums = []
		for i, tad in enumerate(tads):
			for j in range(tad[0], tad[1]):
				tadNums.append(i)
	maxIndex = len(tadNums) - 1

	points_to_add = np.zeros(structure.chrom.getLength(), dtype=bool)	#true if locus should be added
	tracker = Tracker("Identifying loci", structure.chrom.size)

	#find which loci should be added
	with open(path) as listFile:
		for line in listFile:
			line = line.strip().split()
			pos1 = int(line[1])
			pos2 = int(line[4])
			pointNum1 = structure.chrom.getPointNum(pos1)
			pointNum2 = structure.chrom.getPointNum(pos2)
			if pointNum1 and pointNum2:
				tadNum1 = tadNums[min(pointNum1, maxIndex)]
				tadNum2 = tadNums[min(pointNum2, maxIndex)]
				if pointNum1 != pointNum2 and tadNum1 == tadNum2:		#must be in same TAD
					points_to_add[pointNum1] = True
					points_to_add[pointNum2] = True
			tracker.increment()
		listFile.close()

	#create points
	points = np.zeros(structure.chrom.getLength(), dtype=np.object)
	pointNums = np.where(points_to_add == True)[0]
	for i, pointNum in enumerate(pointNums):
		points[pointNum] = Point((0,0,0), pointNum, structure.chrom, i)
	structure.points = points
	
	return structure

def chromFromBed(path):
	"""Initialize ChromParams from intrachromosomal file in BED format"""
	minPos = sys.float_info.max
	maxPos = 0
	print "Scanning {}".format(path)
	with open(path) as infile:
		for i, line in enumerate(infile):
			line = line.strip().split()
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
			if i == 0:
				name = line[0]
				res = (int(line[2]) - pos1)	
		infile.close()
	minPos = int(np.ceil(float(minPos)/res)) * res	#round
	maxPos = int(np.ceil(float(maxPos)/res)) * res
	return ChromParameters(minPos, maxPos, res, name, i)

def basicParamsFromBed(path):
	print "Scanning {}".format(path)
	with open(path) as infile:
		for i, line in enumerate(infile):
			if i == 0:
				line = line.strip().split()
				res = (int(line[2]) - int(line[1]))
		infile.close()
	return i, res

def matFromBed(path, structure):	
	"""Converts BED file to matrix. Only includes loci in structure."""
	if structure is None:
		structure = structureFromBed(path, None, None)

	pointNums = structure.getPointNums()

	numpoints = len(pointNums)
	mat = np.zeros((numpoints, numpoints))	

	maxPointNum = max(pointNums)
	assert maxPointNum - structure.offset < len(structure.points)

	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			loc1 = int(line[1])
			loc2 = int(line[4])
			index1 = structure.getIndex(loc1)
			index2 = structure.getIndex(loc2)
			if index1 and index2:
				if index1 > index2:
					row = index1
					col = index2
				else:
					row = index2
					col = index1
				mat[row, col] += float(line[6])
		infile.close()

	at.makeSymmetric(mat)
	rowsums = np.array([sum(row) for row in mat])
	point_nums = structure.getPointNums()[np.where(rowsums == 0)[0]]
	assert len(np.where(rowsums == 0)[0]) == 0

	return mat

def highToLow(highstructure, resRatio):
	"""Reduces resolution of structure"""
	lowChrom = highstructure.chrom.reduceRes(resRatio)

	low_n = int(np.ceil(len(highstructure.points)/float(resRatio)))

	lowstructure = Structure(np.zeros(low_n, dtype=np.object), [], lowChrom, highstructure.offset/resRatio)

	allPointsToMerge = []
	for i in range(len(lowstructure.points)):
		allPointsToMerge.append([])
	
	for highPoint in highstructure.getPoints():
		pointsToMerge = []
		highNum = highPoint.num
		lowNum = highNum/resRatio
		allPointsToMerge[lowNum - lowstructure.offset].append(highPoint)

	index = lowstructure.offset
	for i, pointsToMerge in enumerate(allPointsToMerge):
		if len(pointsToMerge) > 0:
			lowstructure.points[i] = mergePoints(pointsToMerge, i + lowstructure.offset, lowChrom, index)
			index += 1

	return lowstructure

def mergePoints(pointsToMerge, newPointNum, chrom, index):
	"""Creates new point with average position of pointsToMerge"""
	coords = np.array([point.pos for point in pointsToMerge])
	meanCoord = np.mean(coords, axis=0)
	return Point(meanCoord, newPointNum, chrom, index)

def structure_from_file(path):
	hasMore = True
	with open(path) as infile:
		name = infile.readline().strip()
		res = int(infile.readline().strip())
		minPos = int(infile.readline().strip())
		chrom = ChromParameters(minPos, None, res, name, None)
		structure = Structure([], [], chrom, 0)
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
				structure.points.append(point)
		infile.close()
	structure.points = np.array(structure.points)
	structure.chrom.maxPos = structure.chrom.minPos + structure.chrom.res*num	#max pos is last point num
	return structure
