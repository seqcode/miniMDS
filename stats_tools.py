import numpy as np
import linear_algebra as la

def rmsd(cluster1, cluster2):
	"""Root mean square distance"""
	assert cluster1.chrom.res == cluster2.chrom.res
	res = cluster1.chrom.res
	assert cluster1.chrom.minPos/res == cluster2.chrom.minPos/res	#indexing must be same	

	intersection = [num for num in cluster1.getPointNums() if num in cluster2.getPointNums()]

	dist_sum = 0
	for num in intersection:
		point1 = cluster1.points[num - cluster1.offset]
		point2 = cluster2.points[num - cluster2.offset]
		dist_sum += la.calcDistance(point1.pos, point2.pos)**2
	msd = dist_sum/len(intersection)	#mean square distance
	return msd**(1./2)	#root mean square distance

def movingAverage(signal, size_of_window):
	"""Modified from http://beauty-of-imagination.blogspot.fr/2012/09/fun-with-signal-processing-and.html"""
	window = np.ones(size_of_window)
	return np.roll(np.convolve(window/size_of_window, signal, "valid" ), size_of_window/2)

def smoothWithMovingAverage(signal, size_of_window):
	smoothed = movingAverage(signal, size_of_window)
	signal_size = len(signal)
	remainder = signal[signal_size - size_of_window + 1 : signal_size]	#end of signal, which can't be smoothed
	smoothed_remainder = np.zeros_like(remainder)
	remainder_size = size_of_window - 1
	for i in range(remainder_size):
		smoothed_remainder[i] = movingAverage(remainder[i:remainder_size], remainder_size-i)
	return np.concatenate((smoothed, smoothed_remainder))

def cmds(distMat):
	"""Modified from http://www.nervouscomputer.com/hfs/cmdscale-in-python/"""
	# Number of points                                                                        
	n = len(distMat)

	# Centering matrix                                                                        
	h = np.eye(n) - np.ones((n, n))/n

	# YY^T                                                                                    
	b = -h.dot(distMat**2).dot(h)/2

	# Diagonalize                                                                             
	evals, evecs = np.linalg.eigh(b)

	# Sort by eigenvalue in descending order                                                  
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	evecs = evecs[:,idx]

	return np.array([evecs[:,0]*evals[0]**(1./2), evecs[:,1]*evals[1]**(1./2), evecs[:,2]*evals[2]**(1./2)]).T
