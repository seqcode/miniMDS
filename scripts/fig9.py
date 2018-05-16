from mayavi import mlab
import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot
import linear_algebra as la
import numpy as np

def plot_coords_interactive(coords, res, color=(1,0,0), radius=None, out_path=None):
	if radius is None:
		radius = calculateRadius(coords, res)
	xs = coords[:,0]
	ys = coords[:,1]
	zs = coords[:,2]
	mlab.figure(bgcolor=(1,1,1))
	mlab.plot3d(xs, ys, zs, tube_radius=radius, color=color)
	if out_path is not None:
		mlab.savefig(out_path)	
	mlab.show()

def calculateRadius(coords, res):
	"""Calculate to-scale radius based on Kuhn length and diameter of chromatin"""
	#from Rippe (2001)
	kl = 289	#Kuhn length (nm)
	bpPerKL = 30000.	#base pairs per Kuhn length 
	chromatinDiameter = 30	#diameter of heterochromatin (nm)

	totDist = 0
	count = 0
	n = len(coords)
	for i in range(1, n):
		totDist += la.calcDistance(coords[i-1], coords[i])
		count += 1
	avgDist = totDist/count		#average distance between neighboring loci
	physicalDist = kl * (res/bpPerKL)**(1./2)		#physical distance between neighboring loci (nm)
	conversionFactor = avgDist/physicalDist
	return chromatinDiameter/2 * conversionFactor

mmds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_10kb_mmds_coords.tsv")
cmds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_10kb_cmds_coords.tsv")
minimds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_10kb_minimds_coords.tsv")

mmds_res = mmds_structure.chrom.res
cmds_res = cmds_structure.chrom.res
minimds_res = minimds_structure.chrom.res

assert mmds_res == cmds_res == minimds_res

res = mmds_res

plot.plot_structure_interactive(mmds_structure, out_path="Fig9A.png")
plot.plot_structure_interactive(cmds_structure, out_path="Fig9B.png")
plot.plot_structure_interactive(minimds_structure, out_path="Fig9C.png")
plot_coords_interactive(np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep1_coords.tsv"), res, out_path="Fig9D.png")
