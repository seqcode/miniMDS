import numpy as np
from scipy import stats as st
import sys
sys.path.append("..")
import linear_algebra as la

def distsFromCoords(coords):
	"""Creates distance matrix from 3D coords"""
	n = len(coords)
	distMat = np.zeros((n,n))
	for i in range(n):
		for j in range(i):
			distMat[i,j] = la.calcDistance(coords[i], coords[j])
			if distMat[i,j] == 0:
				print "Error. Duplicate coordinates."
				print coords[i]
				print coords[j]
				sys.exit(0)
	return distMat

def pearson(mat1, mat2):
	"""Root mean square error between two matrices, ignoring zeroes"""
	assert mat1.shape == mat2.shape
	#convert to vectors
	vec1 = mat1.flatten()
	vec2 = mat2.flatten()

	#remove zeroes
	nonzero = [i for i in range(len(vec1)) if vec1[i] != 0 and vec2[i] != 0]
	vec1 = vec1[nonzero]
	vec2 = vec2[nonzero]

	r, p = st.pearsonr(vec1, vec2)
	return r

def parse_time(time_string):
	split = time_string.split("m")
	mins = int(split[0])
	secs = float(split[1].split("s")[0])
	return mins + secs/60

def distMat(cluster):
	"""Creates distance matrix from cluster"""
	points = cluster.getPoints()
	numPoints = len(points)
	mat = np.zeros((numPoints, numPoints))
	for i in range(numPoints):
		for j in range(i):
			mat[i,j] = la.calcDistance(points[i].pos, points[j].pos)
	return mat
