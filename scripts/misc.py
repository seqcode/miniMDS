from mayavi import mlab	#this must be first
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

def plot_coords_interactive(coords, res, color=(1,0,0), radius=None, out_path=None):
	if radius is None:
		radius = calculateRadius(coords, res)
	xs = coords[:,0]
	ys = coords[:,1]
	zs = coords[:,2]
	mlab.figure(bgcolor=(1,1,1))
	mlab.plot3d(xs, ys, zs, tube_radius=radius, color=color)
	if out_path is not None:
		mlab.savefig(out_path)	
	mlab.show()

def calculateRadius(coords, res):
	"""Calculate to-scale radius based on Kuhn length and diameter of chromatin"""
	#from Rippe (2001)
	kl = 289	#Kuhn length (nm)
	bpPerKL = 30000.	#base pairs per Kuhn length 
	chromatinDiameter = 30	#diameter of heterochromatin (nm)

	totDist = 0
	count = 0
	n = len(coords)
	for i in range(1, n):
		totDist += la.calcDistance(coords[i-1], coords[i])
		count += 1
	avgDist = totDist/count		#average distance between neighboring loci
	physicalDist = kl * (res/bpPerKL)**(1./2)		#physical distance between neighboring loci (nm)
	conversionFactor = avgDist/physicalDist
	return chromatinDiameter/2 * conversionFactor
