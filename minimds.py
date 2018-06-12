import sys
import numpy as np
from sklearn import manifold
import argparse
import pymp
import multiprocessing as mp
import data_tools as dt
import array_tools as at
import tad
import linear_algebra as la
import tools

def infer_structure(contactMat, structure, alpha, num_threads, classical=False):
	"""Infers 3D coordinates for one structure"""
	assert len(structure.nonzero_abs_indices()) == len(contactMat)

	at.makeSymmetric(contactMat)
	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)
	at.makeSymmetric(distMat)

	if classical:	#classical MDS
		coords = la.cmds(distMat)
	else:
		coords = manifold.MDS(n_components=3, metric=True, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=num_threads).fit_transform(distMat)

	structure.setCoords(coords)

def fullMDS(path, classical, alpha, num_threads):
	"""MDS without partitioning"""
	structure = dt.structureFromBed(path)
	contactMat = dt.matFromBed(path, structure)
	infer_structure(contactMat, structure, alpha, num_threads, classical)
	return structure
	
def partitionedMDS(path, args):
	"""Partitions structure into substructures and performs MDS"""
	domainSmoothingParameter = args[0]
	minSizeFraction = args[1]
	maxmemory = args[2]
	num_threads = args[3]
	alpha = args[4]
	res_ratio = args[5]
	alpha2 = args[6]

	#create low-res structure
	low_chrom = dt.chromFromBed(path)
	low_chrom.res *= res_ratio
	lowstructure = dt.structureFromBed(path, low_chrom)	#low global structure

	#get TADs
	low_contactMat = dt.matFromBed(path, lowstructure)
	low_tad_indices = tad.getDomains(low_contactMat, lowstructure, domainSmoothingParameter, minSizeFraction)		#low substructures, defined on relative indices not absolute indices
	tad.substructuresFromTads(lowstructure, low_tad_indices)

	#create high-res chrom
	size, res = dt.basicParamsFromBed(path)
	highChrom = dt.ChromParameters(lowstructure.chrom.minPos, lowstructure.chrom.maxPos, res, lowstructure.chrom.name, size)

	highstructure = dt.Structure([], [], highChrom, 0)
	high_substructures = []
	
	low_gen_coords = lowstructure.getGenCoords()
	offset = 0 #initialize
	for td in low_tad_indices:
		start_gen_coord = low_gen_coords[td[0]]
		end_gen_coord = low_gen_coords[td[1]]
		high_substructure = dt.structureFromBed(path, highChrom, start_gen_coord, end_gen_coord, offset)
		high_substructures.append(high_substructure)
		offset += len(high_substructure.points)	#update
		offset -= 1
	
	highstructure.setstructures(high_substructures)
	
	infer_structure(low_contactMat, lowstructure, alpha, num_threads)
	print "Low-resolution MDS complete"

	highSubstructures = pymp.shared.list(highstructure.structures)
	lowSubstructures = pymp.shared.list(lowstructure.structures)

	numSubstructures = len(highstructure.structures)
	num_threads = min((num_threads, mp.cpu_count(), numSubstructures))	#don't exceed number of requested threads, available threads, or structures
	with pymp.Parallel(num_threads) as p:
		for substructurenum in p.range(numSubstructures):
			highSubstructure = highSubstructures[substructurenum]	
			if len(highSubstructure.getPoints()) > 0:	#skip empty
				trueLow = lowSubstructures[substructurenum]

				#perform MDS individually
				structure_contactMat = dt.matFromBed(path, highSubstructure)	#contact matrix for this structure only
				infer_structure(structure_contactMat, highSubstructure, alpha2, num_threads)

				#approximate as low resolution
				inferredLow = dt.highToLow(highSubstructure, res_ratio)

				#rescale
				scaling_factor = la.radius_of_gyration(trueLow)/la.radius_of_gyration(inferredLow)
				for i, point in enumerate(inferredLow.points):
					if point != 0:
						x, y, z = point.pos
						inferredLow.points[i].pos = (x*scaling_factor, y*scaling_factor, z*scaling_factor)

				#recover the transformation for inferred from true low structure
				r, t = la.getTransformation(inferredLow, trueLow)
				t /= scaling_factor

				#transform high structure
				highSubstructure.transform(r, t)
				highSubstructures[substructurenum] = highSubstructure

				print "MDS performed on structure {} of {}".format(substructurenum + 1, numSubstructures)

	highstructure.setstructures(highSubstructures)

	return highstructure

def main():
	parser = argparse.ArgumentParser(description="Reconstruct 3D coordinates from normalized intrachromosomal Hi-C BED files.")
	parser.add_argument("path", help="path to intrachromosomal Hi-C BED file")
	parser.add_argument("--classical", action="store_true", help="use classical MDS (default: metric MDS)")
	parser.add_argument("--full", action="store_true", help="use full MDS (default: partitioned MDS)")
	parser.add_argument("-l", type=int, help="low resolution/high resolution", default=10)
	parser.add_argument("-p", type=float, default=0.1, help="domain size parameter: larger value means fewer structures created (for partitioned MDS only)")
	parser.add_argument("-m", type=float, default=0.05, help="minimum domain size parameter: prevents structures from being too small (for partitioned MDS only)")
	parser.add_argument("-o", help="path to output file")
	parser.add_argument("-r", default=32000000, help="maximum RAM to use (in kb)")
	parser.add_argument("-n", type=int, default=3, help="number of threads")
	parser.add_argument("-a", type=float, default=4, help="alpha factor for converting contact frequencies to physical distances")
	parser.add_argument("-a2", type=float, default=2.5, help="short-range alpha factor for converting contact frequencies to physical distances")
	args = parser.parse_args()

	if args.full:	#not partitioned
		structure = fullMDS(args.path, args.classical, args.a, args.n)

	else:	#partitioned
		params = (args.p, args.m, args.r, args.n, args.a, args.l, args.a2)
		names = ("Domain size parameter", "Minimum domain size", "Maximum memory", "Number of threads", "Alpha", "Resolution ratio", "Short-range, alpha")
		intervals = ((0, 1), (0, 1), (0, None), (0, None), (1, None), (1, None), (1, None))
		if not tools.args_are_valid(params, names, intervals):
			sys.exit(1)

		structure = partitionedMDS(args.path, params)
	
	if args.o:
		structure.write(args.o)
	else:
		prefix = args.path.split(".bed")[0]
		structure.write(prefix + "_structure.tsv")

if __name__ == "__main__":
	main()
