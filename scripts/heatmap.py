import matplotlib
from matplotlib import pyplot as plt
import numpy as np

def threshold(mat, value):
	"""Cuts off values above threshold for ease of visualization of heatmap"""
	n = len(mat)
	thresholded = np.zeros_like(mat)
	for i in range(n):
		for j in range(len(mat[0])):
			thresholded[i,j] = min((mat[i,j], value))
	return thresholded

def createHeatmap(mat, domains, outpath, colors=None):
	# Plot
	fig, ax = plt.subplots()
	plt.pcolormesh(mat, cmap=plt.cm.Reds)

	# turn off the frame
	ax.set_frame_on(False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	#turn off ticks
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	for t in ax.xaxis.get_major_ticks():
		t.tick1On = False
		t.tick2On = False
	for t in ax.yaxis.get_major_ticks():
		t.tick1On = False
		t.tick2On = False

	#plot domain boundaries
	if domains is not None:
		if colors is None:	#default is black
			colors = ["k" for domain in domains]

		for (domain, color) in zip(domains, colors):
			lowerBound = domain[0]
			upperBound = domain[1]
			plt.plot([lowerBound, upperBound], [lowerBound, lowerBound], c=color, lw=1)	#horizontal
			plt.plot([lowerBound, upperBound], [upperBound, upperBound], c=color, lw=1)	#lower horizontal
			plt.plot([upperBound, upperBound], [lowerBound, upperBound], c=color, lw=1)	#vertical
			plt.plot([lowerBound, lowerBound], [lowerBound, upperBound], c=color, lw=1)	#left vertical

	if outpath is not None:
		plt.savefig(outpath)

	else:
		plt.show()

def heatMapFromMat(mat, maxvalue=None, tads=None, outpath=None, colors=None):
	if maxvalue is not None:
		mat = threshold(mat, maxvalue)
	createHeatmap(mat, tads, outpath, colors)
