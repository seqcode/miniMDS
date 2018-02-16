from mayavi import mlab
import numpy as np
import linear_algebra as la
import os
import sys

#from Rippe 2001
kl = 289	#Kuhn length (nm)
bpPerKL = 30000.	#base pairs per Kuhn length 
chromatinDiameter = 30	#diameter of heterochromatin (nm)

default_colors = np.array([[255,0,0], [0,255,238], [255,238,0], [0,102,255], [255,0,170], [255,102,0], [204,255,0], [0,238,255], [0,68,255], [255,0,102], [255,136,0], [0,255,34], [0,204,255], [34,0,255], [255,0,68], [255,170,0], [0,255,136], [0,170,255], [204,0,255], [255,204,0], [0,255,204], [0,136,255], [255,0,238]])/255.
default_colors = [tuple(color) for color in default_colors]	#convert to tuple

def plot_structures_interactive(structures, all_enrichments=None, colors=default_colors, radius=None, cut=False, out_path=None):
	mlab.close(all=True)
	mlab.figure(bgcolor=(1,1,1))
	if radius is None:
		radius = calculateRadius(structures)
	for i, structure in enumerate(structures):
		coords = np.array(structure.getCoords())
		xs = coords[:,0]
		ys = coords[:,1]
		zs = coords[:,2]
		if cut:
			midpoint = np.mean(xs)
			indices = np.where(xs > midpoint)[0]
			xs = xs[indices]
			ys = ys[indices]
			zs = zs[indices]	
		if all_enrichments:
			mlab.plot3d(xs, ys, zs, all_enrichments[i], tube_radius=radius)
		else:
			mlab.plot3d(xs, ys, zs, tube_radius=radius, color=colors[i])
	if out_path:
		mlab.savefig(out_path)		
	mlab.show()

def plot_structure_interactive(structure, enrichments=None, color=(1,0,0), radius=None, out_path=None):
	if radius is None:
		radius = calculateRadius([structure])
	coords = np.array(structure.getCoords())
	xs = coords[:,0]
	ys = coords[:,1]
	zs = coords[:,2]
	mlab.figure(bgcolor=(1,1,1))
	if enrichments:
		mlab.plot3d(xs, ys, zs, tube_radius=radius, color=color)
	else:
		mlab.plot3d(xs, ys, zs, tube_radius=radius, color=color)
	if out_path:
		mlab.savefig(out_path)	
	mlab.show()

def plot_structures_gif(structures, outname, all_enrichments=None, colors=default_colors, radius=None, increment=10):
	if 360%increment != 0:
		print "Error. Increment must be factor of 360."
		sys.exit(0)
	if radius is None:
		radius = calculateRadius(structures)
	mlab.figure(bgcolor=(1,1,1))
	for i in range(0, 360, increment):
		for j, structure in enumerate(structures):
			coords = np.array(structure.getCoords())
			if all_enrichments:
				s = mlab.plot3d(coords[:,0], coords[:,1], coords[:,2], all_enrichments[j], tube_radius=radius)
			else:
				s = mlab.plot3d(coords[:,0], coords[:,1], coords[:,2], tube_radius=radius, color=colors[j])
		mlab.view(i)
		mlab.savefig("{}_{:>03}.png".format(outname, i))

	mlab.close()
	os.system("convert -loop 1 {}_*.png {}.gif".format(outname, outname))
	os.system("rm {}_*.png".format(outname))

def plot_structure_gif(structure, outname, enrichments=None, color=(1,0,0), radius=None, increment=10):
	if 360%increment != 0:
		print "Error. Increment must be factor of 360."
		sys.exit(0)
	if radius is None:
		radius = calculateRadius([structure])
	coords = np.array(structure.getCoords())
	mlab.figure(bgcolor=(1,1,1))
	if enrichments:
		s = mlab.plot3d(coords[:,0], coords[:,1], coords[:,2], enrichments, tube_radius=radius)
	else:
		s = mlab.plot3d(coords[:,0], coords[:,1], coords[:,2], tube_radius=radius, color=color)
	for i in range(0, 360, increment):
		mlab.view(i)
		mlab.savefig("{}_{:>03}.png".format(outname, i))
		
	mlab.close()
	os.system("convert -loop 1 {}_*.png {}.gif".format(outname, outname))
	os.system("rm {}_*.png".format(outname))

def calculateRadius(structures):
	"""Calculate to-scale radius based on Kuhn length and diameter of chromatin"""
	conversionFactors = np.zeros(len(structures))
	for j, structure in enumerate(structures):
		totDist = 0
		count = 0
		coords = structure.getCoords()
		n = len(coords)
		for i in range(1, n):
			totDist += la.calcDistance(coords[i-1], coords[i])
			count += 1
		avgDist = totDist/count		#average distance between neighboring loci
		physicalDist = kl * (structure.chrom.res/bpPerKL)**(1./2)		#physical distance between neighboring loci (nm)
		conversionFactors[j] = avgDist/physicalDist
	conversionFactor = np.mean(conversionFactors)
	return chromatinDiameter/2 * conversionFactor
