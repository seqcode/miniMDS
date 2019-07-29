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

	def getAbsoluteIndex(self, genCoord):
		"""Converts genomic coordinate into absolute index. Absolute indexing includes empty (zero) points."""
		if genCoord < self.minPos or genCoord > self.maxPos + self.res:
			return None
		else:
			return int((genCoord - self.minPos)/self.res) 

	def getGenCoord(self, abs_index):
		"""Converts absolute index into genomic coordinate"""
		return self.minPos + self.res * abs_index

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
		if len(structures) == 0 or structures is None:
			self.structures = []
		else:
			self.setstructures(structures)
		self.chrom = chrom	#chromosome parameters
		self.offset = offset	#absolute indexing offset (for substructures only)

	def getCoords(self):
		return [point.pos for point in self.getPoints()]

	def setCoords(self, coords):
		for coord, abs_index in zip(coords, self.nonzero_abs_indices()):
			self.points[abs_index - self.offset].pos = coord

	def nonzero_abs_indices(self):
		"""Absolute indices for all non-zero points."""
		return np.array([point.absolute_index for point in self.getPoints()])

	def getPoints(self):
		"""All non-zero points"""
		return self.points[np.where(self.points != 0)[0]]

	def subsamplePoints(self, start_abs_index, end_abs_index):
		"""Set structure's points to only include start_abs_index through end_abs_index"""
		points = self.points[start_abs_index:end_abs_index+1]
		self.chrom.maxPos = self.chrom.getGenCoord(end_abs_index)
		self.chrom.minPos = self.chrom.getGenCoord(start_abs_index)
		#re-index
		for abs_index in np.where(points != 0)[0]:
			points[abs_index].absolute_index = abs_index
		self.points = points
		self.set_rel_indices()

	def getGenCoords(self):
		"""Non-zero genomic coordinates of structure"""
		return [self.chrom.getGenCoord(abs_index) for abs_index in self.nonzero_abs_indices()]

	def get_rel_index(self, genCoord):
		"""Converts genomic coordinate into relative index."""
		abs_index = self.chrom.getAbsoluteIndex(genCoord)
		if abs_index is None:
			return None
		else:
			abs_index -= self.offset
			if abs_index >= 0 and abs_index < len(self.points):
				point = self.points[abs_index]
				if point == 0:
					return None
				else:
					return point.relative_index
			else:
				return None
	
	def setstructures(self, structures):
		self.structures = structures
		self.points = np.zeros(max([max(structure.nonzero_abs_indices()) for structure in structures]) + 1, dtype=np.object)	#reset
		for structure in self.structures:
			for point in structure.points:
				if point != 0:
					self.points[point.absolute_index] = point

	def createSubstructure(self, points, offset):
		"""Creates substructure containing points"""
		substructure = Structure(points, [], self.chrom, offset)
		substructure.set_rel_indices()
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
		for i, abs_index in enumerate(self.nonzero_abs_indices()):
			self.points[abs_index - self.offset].pos = a_transformed[i]

	def write(self, outpath):
		with open(outpath, "w") as out:
			out.write(self.chrom.name + "\n")
			out.write(str(self.chrom.res) + "\n")
			out.write(str(self.chrom.minPos) + "\n")
			abs_index = self.offset
			for point in self.points:
				if point == 0:
					out.write("\t".join((str(abs_index), "nan", "nan", "nan")) + "\n")
				else:
					out.write("\t".join((str(abs_index), str(point.pos[0]), str(point.pos[1]), str(point.pos[2]))) + "\n")
				abs_index += 1
		out.close()

	def set_rel_indices(self):
		"""Relative indexing is index relative to non-zero points only"""
		for i, abs_index in enumerate(self.nonzero_abs_indices()):
			self.points[abs_index - self.offset].relative_index = i

	def rescale(self):
		"""Rescale radius of gyration of structure to 1"""
		rg = la.radius_of_gyration(self)
		for i, point in enumerate(self.points):
			if point != 0:
				x, y, z = point.pos
				self.points[i].pos = (x/rg, y/rg, z/rg)

class Point(object):
	"""Point in 3-D space"""
	def __init__(self, pos, chrom, absolute_index, relative_index):
		self.pos = pos	#3D coordinates
		self.chrom = chrom	#chromosome parameters
		self.absolute_index = absolute_index	#index relative to all points in structure (including null/zero points)
		self.relative_index = relative_index	#index relative to only non-zero points

def structureFromBed(path, chrom=None, start=None, end=None, offset=0):
	"""Initializes structure from intrachromosomal BED file."""
	if chrom is None:
		chrom = chromFromBed(path)

	if start is None:
		start = chrom.minPos

	if end is None:
		end = chrom.maxPos

	structure = Structure([], [], chrom, offset)
	
	structure.points = np.zeros(int((end - start)/chrom.res) + 1, dtype=object)	#true if locus should be added
	tracker = Tracker("Identifying loci", structure.chrom.size)

	#add loci
	with open(path) as listFile:
		for line in listFile:
			line = line.strip().split()
			pos1 = int(line[1])
			pos2 = int(line[4])
			if pos1 >= start and pos1 <= end and pos2 >= start and pos2 <= end:
				abs_index1 = structure.chrom.getAbsoluteIndex(pos1)
				abs_index2 = structure.chrom.getAbsoluteIndex(pos2)
				if abs_index1 != abs_index2:	#non-self-interacting
					structure.points[int((pos1 - start)/chrom.res)] = Point((0,0,0), structure.chrom, abs_index1, 0)
					structure.points[int((pos2 - start)/chrom.res)] = Point((0,0,0), structure.chrom, abs_index2, 0)
			tracker.increment()
		listFile.close()

	structure.set_rel_indices()
	
	return structure

def chromFromBed(path):
	"""Initialize ChromParams from intrachromosomal file in BED format"""
	minPos = sys.float_info.max
	maxPos = 0
	print("Scanning {}".format(path))
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
	minPos = int(np.floor(float(minPos)/res)) * res	#round
	maxPos = int(np.ceil(float(maxPos)/res)) * res
	return ChromParameters(minPos, maxPos, res, name, i)

def basicParamsFromBed(path):
	print("Scanning {}".format(path))
	with open(path) as infile:
		for i, line in enumerate(infile):
			if i == 0:
				line = line.strip().split()
				res = (int(line[2]) - int(line[1]))
		infile.close()
	return i, res

def matFromBed(path, structure=None):	
	"""Converts BED file to matrix. Only includes loci in structure."""
	if structure is None:
		structure = structureFromBed(path, None, None)

	abs_indices = structure.nonzero_abs_indices()

	numpoints = len(abs_indices)
	mat = np.zeros((numpoints, numpoints))	

	assert max(abs_indices) - structure.offset < len(structure.points)

	with open(path) as infile:
		for line in infile:
			line = line.strip().split()
			loc1 = int(line[1])
			loc2 = int(line[4])
			index1 = structure.get_rel_index(loc1)
			index2 = structure.get_rel_index(loc2)
			if index1 is not None and index2 is not None:
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
	if len(np.where(rowsums == 0)[0]) != 0:
		print("Error: empty rows")
		print(np.array(structure.getGenCoords())[np.where(rowsums == 0)[0]])
		sys.exit(1)
	assert len(np.where(rowsums == 0)[0]) == 0

	return mat

def highToLow(highstructure, resRatio):
	"""Reduces resolution of structure"""
	lowChrom = highstructure.chrom.reduceRes(resRatio)

	low_n = int(len(highstructure.points)/resRatio + 1)

	lowstructure = Structure(np.zeros(low_n, dtype=np.object), [], lowChrom, int(highstructure.offset/resRatio))

	allPointsToMerge = [[] for i in range(low_n)]
	
	for highPoint in highstructure.getPoints():
		pointsToMerge = []
		high_abs_index = highPoint.absolute_index - highstructure.offset
		low_abs_index = int(high_abs_index/resRatio)
		allPointsToMerge[low_abs_index].append(highPoint)

	index = lowstructure.offset
	for i, pointsToMerge in enumerate(allPointsToMerge):
		if len(pointsToMerge) > 0:
			meanCoord = np.mean(np.array([point.pos for point in pointsToMerge]), axis=0)
			lowstructure.points[i] = Point(meanCoord, lowChrom, i + lowstructure.offset, index)
			index += 1

	return lowstructure

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
					point = Point((x,y,z), chrom, num, index)
					index += 1
				structure.points.append(point)
		infile.close()
	structure.points = np.array(structure.points)
	structure.chrom.maxPos = structure.chrom.minPos + structure.chrom.res*num	#max pos is last point num
	return structure

def make_compatible(structures):
	"""Enforce that points be shared by all structures"""
	gen_coord_dict = {}
	for i, structure in enumerate(structures):
		for gen_coord in structure.getGenCoords():
			if gen_coord in gen_coord_dict:
				gen_coord_dict[gen_coord] += 1
			else:
				gen_coord_dict[gen_coord] = 1
	
	consensus = []
	n = len(structures)
	for gen_coord in gen_coord_dict.keys():
		if gen_coord_dict[gen_coord] == n:
			consensus.append(gen_coord)

	consensus = np.sort(consensus)
	
	for structure in structures:
		new_chrom = ChromParameters(consensus[0], consensus[-1] + structure.chrom.res, structure.chrom.res, structure.chrom.name, structure.chrom.size)
		new_points = np.zeros(new_chrom.getLength(), dtype=object)
		for i, gen_coord in enumerate(consensus):
			old_abs_index = structure.chrom.getAbsoluteIndex(gen_coord)
			new_abs_index = new_chrom.getAbsoluteIndex(gen_coord)
			pos = structure.points[old_abs_index].pos
			new_points[new_abs_index - structure.offset] = Point(pos, new_chrom, new_abs_index, i)
		structure.points = new_points
		structure.chrom = new_chrom

def normalized_dist_mat(path, structure):
	"""Standard processing for creating distance matrix"""
	contacts = matFromBed(path, structure)
	dists = at.contactToDist(contacts, 4)
	at.makeSymmetric(dists)
	return dists/np.mean(dists)	#normalize
