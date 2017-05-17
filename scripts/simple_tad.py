import sys
sys.path.append("..")
import stats_tools as st
import numpy as np

def calcScore(locNum, contactMat, numLoc):
	"""Calculate directionality index for locus. See Dixon 2012 supplemental."""
	a = 0		#initialize
	b = 0
	aCount = 0
	bCount = 0
	avg_a = 0
	avg_b = 0
	for i in range(locNum-numLoc, locNum):	#upstream
		a += contactMat[locNum][i]
		aCount += 1
	if aCount != 0:
		avg_a=a/aCount
	for i in range(locNum+1, locNum+numLoc):	#downstream
		b += contactMat[i][locNum]
		bCount += 1
	if bCount != 0:
		avg_b = b/bCount
	if avg_a + avg_b != 0 and avg_a != avg_b:
		e = (avg_a + avg_b)/2
		index = (avg_b - avg_a)/abs(avg_b - avg_a)*((avg_a - e)**2/e + (avg_b - e)**2/e)
	else:
		index = 0
	return index

def allScores(contactMat, maxNumLoc):
	"""Calculate all directionality indices for chromosome"""
	dirIndices=[]	
	totNumLoc = len(contactMat)
	for i in range(totNumLoc):
		numLoc = min((maxNumLoc, totNumLoc - i, i))	#avoid going out of range of contact matrix 
		dirIndices.append([calcScore(i, contactMat, numLoc)][0])
	dirIndices = np.array(dirIndices)	
	return dirIndices

def domainsFromScores(indices, minSizeFraction):
	"""Identify domain starts and ends from directionality indices"""
	numLoc = len(indices)
	minNumLoc = minSizeFraction*numLoc
	starts = []
	ends = []
	prevIndex = np.nan	#initialization 
	currNum = 0
	nextstart = 0
	currend = 0
	for index in indices:
		if index > 0:	#downstream bias
			if prevIndex < 0 or prevIndex is np.nan:
				currstart = nextstart	#start of current domain
				currend = currNum	#end of current domain
				nextstart = currNum		#start of next domain
				D = currend + 8
				if 8==D or currend - currstart > minNumLoc:
					starts.append(nextstart)	
				if prevIndex < 0:	#previous is upstream
					if currend - currstart > minNumLoc:
						ends.append(currend)	
		prevIndex = index
		currNum+=1
	domains = []
	i = 0	#index of starts
	j = 0	#index of ends
	if starts[i] > 0:	#if first start isn't 0
		domains.append((0, ends[j]))
		j+=1
	numstarts = len(starts)
	numends = len(ends)
	while i<numstarts and j<numends:
		domains.append((starts[i], ends[j]))
		i+=1
		j+=1
	if ends[numends-1] < numLoc:	#if last end isn't end of indices
		domains.append((starts[numstarts-1], len(indices)))
	return domains

def getDomains(contactMat, smoothingFactor, minSizeFraction):
	scores = allScores(contactMat, 50)
	smoothed = st.smoothWithMovingAverage(scores, smoothingFactor)
	return domainsFromScores(smoothed, minSizeFraction)
