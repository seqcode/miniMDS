import numpy as np

def calcScore(pointNum, points, contactMat, numPoints):
	"""Calculates directionality score for locus. See Dixon 2012 supplemental. Positive=downstream. Negative=upstream."""
	a = 0	#initialize
	b = 0
	aCount = 0
	bCount = 0
	avg_a = 0
	avg_b = 0
	currIndex = points[pointNum].index

	#upstream
	upstreamIndexFound = False
	numPointsUpstream = numPoints	#numPoints is max separation from curr point to include when calculating score
	while not upstreamIndexFound and numPointsUpstream > 0:
		if points[pointNum - numPointsUpstream] != 0:
			upstreamIndexFound = True
			minIndex = points[pointNum - numPointsUpstream].index	
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
		if points[pointNum + numPointsDownstream] != 0:
			downstreamIndexFound = True
			maxIndex = points[pointNum + numPointsDownstream].index	
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
	for pointNum in structure.getPointNums():
		i = structure.points[pointNum - structure.offset].index
		numPoints = min((maxNumPoints, totNumLoc - 1 - i, i))	#avoid going out of range of contact matrix 
		scores.append(calcScore(pointNum, structure.points, contactMat, numPoints))
	return scores

def domainsFromScores(scores, minSizeFraction):
	minNumLoc = minSizeFraction*len(scores)
	start = 0	#initialize
	prevScore = 0
	domains = []
	for i in range(len(scores)):
		score = scores[i]
		if i == len(scores) - 1:
			end = i + 1
			domains.append([start,end])
			start = i + 1
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

def smoothWithMovingAverage(signal, size_of_window):
	"""Modified from http://beauty-of-imagination.blogspot.fr/2012/09/fun-with-signal-processing-and.html"""
	window = np.ones(size_of_window)
	smoothed = np.roll(np.convolve(window/size_of_window, signal, "valid" ), size_of_window/2)
	signal_size = len(signal)
	remainder = signal[signal_size - size_of_window + 1 : signal_size]	#end of signal, which can't be smoothed
	smoothed_remainder = np.zeros_like(remainder)
	remainder_size = size_of_window - 1
	for i in range(remainder_size):
		smoothed_remainder[i] = movingAverage(remainder[i:remainder_size], remainder_size-i)
	return np.concatenate((smoothed, smoothed_remainder))

def substructuresFromTads(high_structure, low_structure, low_tads):
	"""Create compatible substructures from TADs in high-res structure and low-res structure"""
	res_ratio = low_structure.chrom.res/high_structure.chrom.res
	high_tads = low_tads * res_ratio
        high_offset = 0
	low_offset = 0
        for high_tad, low_tad in zip(high_tads, low_tads):
		high_start = high_tad[0]
		high_end = high_tad[1]
		low_start = low_tad[0]
		low_end = low_tad[1]
                high_points = high_structure.points[(high_start-high_structure.offset):(high_end-high_structure.offset)]
		low_points = low_structure.points[(low_start-low_structure.offset):(low_end-low_structure.offset)]
		high_nums = [high_point.num for high_point in high_points if high_point != 0]
		inferred_low_nums = np.array(high_nums)/res_ratio
		true_low_nums = [low_point.num for low_point in low_points if low_point != 0]
		intersection = [num for num in true_low_nums if num in inferred_low_nums]
		if len(intersection) > 0:	#compatible
                	high_structure.createSubstructure(high_points, high_offset)
			low_structure.createSubstructure(low_points, low_offset)
                high_offset = high_end
		low_offset = low_end
