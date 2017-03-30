import numpy as np
from sklearn import manifold
import sys
sys.path.append("../minimds")
import array_tools as at
import linear_algebra as la

def distToContact(distMat, alpha=-1./3):
	"""Convert distance matrix to contact matrix. Matrix must be symmetric."""
	contactMat = np.zeros_like(distMat)
	numRows = len(contactMat)
	for i in range(numRows):
		for j in range(i+1):
			if distMat[i,j] != 0:
				contactMat[i,j] = distMat[i,j]**(1./alpha)		#see Varoquaux, et al (2014)
	return contactMat

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

def mds(distMat, metric):
	"""Performs MDS on distance matrix"""
	distMat = removeInfinite(distMat)
	mds = manifold.MDS(n_components=3, metric=metric, random_state=np.random.RandomState(seed=3), verbose=0, dissimilarity="precomputed", n_jobs=-1)
	return mds.fit(distMat).embedding_

def removeInfinite(mat):
	"""Replaces infinite values in matrix with zeroes"""
	n = len(mat)
	copy = np.copy(mat)
	for i in range(n):
		for j in range(i+1):
			if not np.isfinite(copy[i,j]):
				 copy[i,j] = 0
	at.makeSymmetric(copy)
	return copy

def pearson(true, embedded):
	"""Root mean square error between two matrices, ignoring zeroes"""
	assert true.shape == embedded.shape
	#convert to vectors
	true = true.flatten()
	embedded = embedded.flatten()

	#remove zeroes
	indices = np.where(true != 0)[0]
	true = true[indices]
	embedded = embedded[indices]

	r, p = st.pearsonr(true, embedded)
	return r
