import numpy as np

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
