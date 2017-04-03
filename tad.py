import stats_tools as st
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

def allScores(contactMat, cluster, maxNumPoints):
	"""Calculates all directionality scores for chromosome"""
	scores = []	
	totNumLoc = len(contactMat)
	for pointNum in cluster.getPointNums():
		i = cluster.points[pointNum - cluster.offset].index
		numPoints = min((maxNumPoints, totNumLoc - 1 - i, i))	#avoid going out of range of contact matrix 
		scores.append(calcScore(pointNum, cluster.points, contactMat, numPoints))
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

def getDomains(contactMat, cluster, sizeParameter, minSizeFraction):
	"""Identify TADs in contact matrix"""
	scores = allScores(contactMat, cluster, 50)	#50 is from Dixon 2012 supplemental
	smoothingFactor = max((int(len(contactMat)*sizeParameter), 1))	#must be >= 1
	smoothed = st.smoothWithMovingAverage(scores, smoothingFactor)
	return domainsFromScores(smoothed, minSizeFraction)

def subclustersFromTads(high_cluster, low_cluster, low_tads):
	"""Create compatible subclusters from TADs in high-res cluster and low-res cluster"""
	res_ratio = low_cluster.chrom.res/high_cluster.chrom.res
	high_tads = low_tads * res_ratio
        high_offset = 0
	low_offset = 0
        for high_tad, low_tad in zip(high_tads, low_tads):
		high_start = high_tad[0]
		high_end = high_tad[1]
		low_start = low_tad[0]
		low_end = low_tad[1]
                high_points = high_cluster.points[(high_start-high_cluster.offset):(high_end-high_cluster.offset)]
		low_points = low_cluster.points[(low_start-low_cluster.offset):(low_end-low_cluster.offset)]
		high_nums = [high_point.num for high_point in high_points if high_point != 0]
		inferred_low_nums = np.array(high_nums)/res_ratio
		true_low_nums = [low_point.num for low_point in low_points if low_point != 0]
		intersection = [num for num in true_low_nums if num in inferred_low_nums]
		if len(intersection) > 0:	#compatible
                	high_cluster.createSubcluster(high_points, high_offset)
			low_cluster.createSubcluster(low_points, low_offset)
                high_offset = high_end
		low_offset = low_end
