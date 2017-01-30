import numpy as np
import stats_tools as st
from tools import Tracker

def contactToDist(contactMat, alpha=-1./3):
	"""Convert contact matrix to distance matrix. Matrix must be symmetric."""
	distMat = np.zeros_like(contactMat)
	numRows = len(contactMat)
	for i in range(numRows):
		for j in range(i+1):
			if contactMat[i,j] != 0:
				distMat[i,j] = contactMat[i,j]**alpha		#see Varoquaux, et al (2014)
	return distMat

def makeSymmetric(mat):
	"""Make below-diagonal matrix symmetric, in place"""
	for row in range(len(mat)):
		for col in range(row):
			mat[col][row] = mat[row][col]

#def sp_interpolate(mat):
#	"""Shortest-path interpolation. See Lesne, et al (2014). Must be intrachromosomal."""
#	mat[np.where(mat == 0)] = np.inf
#	numrows = len(mat)	
#	numcols = len(mat[0])
#	assert numrows == numcols
#	tracker = Tracker("Shortest-path interpolation", numrows)
#	for k in range(numrows):
#		row = mat[k]
#		col = np.array([mat[:,k]]).T	
#		i2k = np.tile(col, [1, numrows])
#		k2j = np.tile(row, [numrows, 1])
#		mat = minArray(mat, i2k+k2j)
#		#tracker.increment()
#	return mat

def minArray(a, b):
	"""Equivalent to MATLAB min() method"""
	assert a.shape == b.shape
	c = np.zeros_like(a)
	for i in range(len(a)):
		for j in range(len(a[0])):
			c[i,j] = min((a[i,j], b[i,j]))
	return c
