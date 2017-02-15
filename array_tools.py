import numpy as np

def contactToDist(contactMat, alpha=-1./3):
	"""Convert contact matrix to distance matrix."""
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
