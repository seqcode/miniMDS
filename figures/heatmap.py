import sys
sys.path.append("..")
from matplotlib import pyplot as plt
import numpy as np
import array_tools as at

def threshold(mat, value):
	"""Cuts off values above threshold for ease of visualization of heatmap"""
	n = len(mat)
	for i in range(n):
		for j in range(n):
			if mat[i,j] > value:
				mat[i,j] = value

def createHeatmap(mat, domains, outpath, colors=None):
	# Plot
	fig, ax = plt.subplots()
	heatmap = ax.pcolor(mat, cmap=plt.cm.Reds)

	# Format
	fig = plt.gcf()

	# turn off the frame
	ax.set_frame_on(False)

	# put the major ticks at the middle of each cell
	ax.set_yticks(np.arange(mat.shape[0]) + 0.5, minor=False)
	ax.set_xticks(np.arange(mat.shape[1]) + 0.5, minor=False)
	ax.set_xticklabels([])
	ax.set_yticklabels([])

	# want a more natural, table-like display
	ax.invert_yaxis()
	ax.xaxis.tick_top()

	# rotate the ticks
	plt.xticks(rotation=90)

	ax.grid(False)

	# Turn off all the ticks
	ax = plt.gca()

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

		for (domain, color) in zip((domains, colors)):
			lowerBound = domain[0]
			upperBound = domain[1]
			plt.plot([lowerBound, upperBound], [lowerBound, lowerBound], c=color, lw=5)	#horizontal
			plt.plot([lowerBound, upperBound], [upperBound, upperBound], c=color, lw=5)	#lower horizontal
			plt.plot([upperBound, upperBound], [lowerBound, upperBound], c=color, lw=5)	#vertical
			plt.plot([lowerBound, lowerBound], [lowerBound, upperBound], c=color, lw=5)	#left vertical
	
	if outpath is not None:
		plt.savefig(outpath)
	plt.show()

def heatMapFromMat(mat, maxvalue, tads, outpath):
	at.makeSymmetric(mat)
	if maxvalue is not None:
		threshold(mat, maxvalue)
	createHeatmap(mat, tads, outpath)
