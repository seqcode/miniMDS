import numpy as np

def calcScore(abs_index, points, contactMat, numPoints):
	"""Calculates directionality score for locus. See Dixon 2012 supplemental. Positive=downstream. Negative=upstream."""
	a = 0	#initialize
	b = 0
	aCount = 0
	bCount = 0
	avg_a = 0
	avg_b = 0
	currIndex = points[abs_index].relative_index

	#upstream
	upstreamIndexFound = False
	numPointsUpstream = numPoints	#numPoints is max separation from curr point to include when calculating score
	while not upstreamIndexFound and numPointsUpstream > 0:
		if points[abs_index - numPointsUpstream] != 0:
			upstreamIndexFound = True
			minIndex = points[abs_index - numPointsUpstream].relative_index	
		else:
			numPointsUpstream -= 1

	if upstreamIndexFound:
		for i in range(minIndex, currIndex):	
			a += contactMat[currIndex][i]
			aCount += 1
		if aCount != 0:
			avg_a = a/aCount

	#downstream
	downstreamIndexFound = False
	numPointsDownstream = numPoints	#numPoints is max separation from curr point to include when calculating score
	while not downstreamIndexFound and numPointsDownstream > 0:
		if points[abs_index + numPointsDownstream] != 0:
			downstreamIndexFound = True
			maxIndex = points[abs_index + numPointsDownstream].relative_index	
		else:
			numPointsDownstream -= 1

	if downstreamIndexFound:
		for i in range(currIndex + 1, maxIndex):	
			b += contactMat[i][currIndex]
			bCount += 1
		if bCount != 0:
			avg_b = b/bCount

	if aCount != 0 and bCount != 0 and avg_a != avg_b:
		e = (avg_a + avg_b)/2
		score = (avg_b-avg_a)/abs(avg_b-avg_a)*((avg_a-e)**2/e+(avg_b-e)**2/e)
	else:
		score = 0

	return score 

def allScores(contactMat, structure, maxNumPoints):
	"""Calculates all directionality scores for chromosome"""
	scores = []	
	totNumLoc = len(contactMat)
	for abs_index in structure.nonzero_abs_indices():
		i = structure.points[abs_index - structure.offset].relative_index
		numPoints = min((maxNumPoints, totNumLoc - 1 - i, i))	#avoid going out of range of contact matrix 
		scores.append(calcScore(abs_index, structure.points, contactMat, numPoints))
	return scores

def domainsFromScores(scores, minSizeFraction):
	minNumLoc = minSizeFraction*len(scores)
	start = 0	#initialize
	prevScore = 0
	domains = []
	for i in range(len(scores)):
		score = scores[i]
		if i == len(scores) - 1:	#at end of scores, close remaining TAD
			end = i
			domains.append([start,end])
		elif score > 0 and prevScore < 0 and i-start >= minNumLoc:	#current is downstream, previous was upstream
			end = i
			domains.append([start,end])
			start = i
		prevScore = score
	return np.array(domains)

def getDomains(contactMat, structure, sizeParameter, minSizeFraction):
	"""Identify TADs in contact matrix"""
	scores = allScores(contactMat, structure, 50)	#50 is from Dixon 2012 supplemental
	smoothingFactor = max((int(len(contactMat)*sizeParameter), 1))	#must be >= 1
	smoothed = smoothWithMovingAverage(scores, smoothingFactor)
	return domainsFromScores(smoothed, minSizeFraction)

def movingAverage(signal, size_of_window):
	"""Modified from http://beauty-of-imagination.blogspot.fr/2012/09/fun-with-signal-processing-and.html"""
	window = np.ones(size_of_window)
	return np.roll(np.convolve(window/size_of_window, signal, "valid"), int(size_of_window/2))

def smoothWithMovingAverage(signal, size_of_window):
	smoothed = movingAverage(signal, size_of_window)
	signal_size = len(signal)
	remainder = signal[signal_size - size_of_window + 1 : signal_size]	#end of signal, which can't be smoothed
	smoothed_remainder = np.zeros_like(remainder)
	remainder_size = size_of_window - 1
	for i in range(remainder_size):
		smoothed_remainder[i] = movingAverage(remainder[i:remainder_size], remainder_size-i)
	return np.concatenate((smoothed, smoothed_remainder))

def substructuresFromTads(structure, tads):
	abs_indices = structure.nonzero_abs_indices()
	offset = 0	#initialize
	for td in tads:
		start = abs_indices[td[0]]	#convert from relative index to absolute index
		end = abs_indices[td[1]]
		points = structure.points[start-structure.offset:end-structure.offset]
		structure.createSubstructure(points, offset)
		offset = end	#update
