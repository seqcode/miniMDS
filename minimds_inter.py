import numpy as np
import data_tools as dt
import linear_algebra as la
import sys
from sklearn import manifold
import tools
import argparse
import minimds as mm
import array_tools as at

def infer_structures(contactMat, structures, offsets, alpha, num_threads, classical=False):
	"""Infers 3D coordinates for multiple structures with same contact matrix"""
	assert sum([len(structure.nonzero_abs_indices()) for structure in structures]) == len(contactMat)

	#at.makeSymmetric(contactMat)
	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)
	at.makeSymmetric(distMat)
	#rowsums = np.array([sum(row) for row in distMat])
	#assert len(np.where(rowsums == 0)[0]) == 0 

	if classical:	#classical MDS
		coords = la.cmds(distMat)
	else:
		coords = manifold.MDS(n_components=3, metric=True, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=num_threads).fit_transform(distMat)

	for offset, structure in zip(offsets, structures):
		structure.setCoords(coords[offset:offset+len(structure.getPoints())])

def get_inter_mat(prefix, inter_res_string, intra_res_string, structures, offsets):
	names = [structure.chrom.name for structure in structures]
	n = len(names)
	for i in range(n):
		if names[i].startswith("chr"):
			names[i] = names[i][3:len(names[i])]	#remove "chr"

	#fill matrix
	total_len = sum([len(structure.getPoints()) for structure in structures])
	mat = np.zeros((total_len, total_len))
	for i in range(n):
		for j in range(i+1):
			if i == j:
				path = "{}_{}_{}.bed".format(prefix, names[i], intra_res_string)
			else:
				path = "{}_{}_{}_{}.bed".format(prefix, names[j], names[i], inter_res_string)
			print("Reading {}".format(path))
			with open(path) as bed:
				for line in bed:
					line = line.strip().split()
					loc1 = int(line[4])
					loc2 = int(line[1])
					index1 = structures[i].get_rel_index(loc1)
					index2 = structures[j].get_rel_index(loc2)
					if index1 is not None and index2 is not None:
						row = index1 + offsets[i]
						col = index2 + offsets[j]
						count = float(line[6])
						mat[row, col] += count
						mat[col, row] += count
			bed.close()
	return mat

def interMDS(names, prefix, inter_res, intra_res, partitioned, args):
	inter_res_string = tools.get_res_string(inter_res)
	intra_res_string = tools.get_res_string(intra_res)

	#get low-res structures from intra files
	low_structures = []
	for name in names:
		path = "{}_{}_{}.bed".format(prefix, name, intra_res_string)
		chrom = dt.chromFromBed(path)
		#reduce res
		chrom.res = inter_res
		chrom.minPos = int(np.floor(float(chrom.minPos)/chrom.res)) * chrom.res	#round
		chrom.maxPos = int(np.ceil(float(chrom.maxPos)/chrom.res)) * chrom.res
		low_structures.append(dt.structureFromBed(path, chrom=chrom))

	#for correct indexing
	n = len(names)
	offsets = np.zeros(n, dtype=int)
	for i in range(1, n):
		offsets[i] = offsets[i-1] + len(low_structures[i-1].getPoints())

	inter_mat = get_inter_mat(prefix, inter_res_string, intra_res_string, low_structures, offsets)

	#perform MDS at low resolution on all chroms
	infer_structures(inter_mat, low_structures, offsets, args[3], args[4])

	#perform MDS at high resolution on each chrom
	high_structures = []
	inferred_low_structures = []
	ts = []
	for true_low, name in zip(low_structures, names):
		path = "{}_{}_{}.bed".format(prefix, name, intra_res_string)
		if partitioned:
			high_structure = mm.partitionedMDS(path, args)
		else:
			high_structure = mm.fullMDS(path, False, args[4], args[3], args[7])
		high_structures.append(high_structure)
		inferred_low = dt.highToLow(high_structure, true_low.chrom.res/high_structure.chrom.res)
		inferred_low_structures.append(inferred_low)

		#rescale
		rescaling_factor = la.radius_of_gyration(true_low)/la.radius_of_gyration(inferred_low)
		rescaled_coords = [rescaling_factor * coord for coord in inferred_low.getCoords()]
		for i, point in enumerate(inferred_low.getPoints()):
			point.pos = rescaled_coords[i]

		r, t = la.getTransformation(inferred_low, true_low)
		high_structure.transform(r, None)	#do not translate now (need to rescale)
		ts.append(t)	

	#translate (with rescaling)
	low_rgs = np.array([la.radius_of_gyration(structure) for structure in low_structures])
	high_rgs = np.array([la.radius_of_gyration(structure) for structure in high_structures])
	scaling_factor = np.mean(high_rgs/low_rgs)
	for high_structure, t in zip(high_structures, ts):
		high_structure.transform(None, scaling_factor*t)	#rescale translation

	return high_structures

def main():
	parser = argparse.ArgumentParser(description="Reconstruct 3D coordinates from normalized interchromosomal Hi-C BED files.")
	parser.add_argument("prefix", help="prefix of Hi-C BED files")
	parser.add_argument("inter_res", type=int, help="resolution of interchromosomal BED files (bp)")
	parser.add_argument("intra_res", type=int, help="resolution of intrachromosomal BED files (bp)")
	parser.add_argument("--partitioned", action="store_true", help="use partitioned MDS (default: full MDS)")
	parser.add_argument("-c", action="append", default=[], help="names of chromosomes to use, e.g. 1 (default: all human chromosomes other than Y)")
	parser.add_argument("-C", type=int, help="number of autosomes")
	parser.add_argument("-l", type=int, help="low resolution/high resolution", default=10)
	parser.add_argument("-p", type=float, default=0.1, help="domain size parameter: larger value means fewer structures created (for partitioned MDS only)")
	parser.add_argument("-m", type=float, default=0.05, help="minimum domain size parameter: prevents structures from being too small (for partitioned MDS only)")
	parser.add_argument("-o", help="prefix of output file")
	parser.add_argument("-r", default=32000000, help="maximum RAM to use (in kb)")
	parser.add_argument("-n", type=int, default=3, help="Number of threads")
	parser.add_argument("-a", type=float, default=4, help="alpha factor for converting contact frequencies to physical distances")
	parser.add_argument("-a2", type=float, default=2.5, help="short-range alpha factor for converting contact frequencies to physical distances")
	parser.add_argument("-w", type=float, default=0.05, help="weight of distance decay prior")
	args = parser.parse_args()

	if len(args.c) == 0:
		if args.C:
			chrom_names = range(1, args.C+1)
		else:
			chrom_names = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
	else:
		chrom_names = args.c
	

	params = (args.p, args.m, args.r, args.n, args.a, args.l, args.a2, args.w)
	names = ("Domain size parameter", "Minimum domain size", "Maximum memory", "Number of threads", "Alpha", "Resolution ratio", "Short-range alpha", "Weight")
	intervals = ((0, 1), (0, 1), (0, None), (0, None), (1, None), (1, None), (1, None), (None, 1))
	if not tools.args_are_valid(params, names, intervals):
		sys.exit(1)

	structures = interMDS(chrom_names, args.prefix, args.inter_res, args.intra_res, args.partitioned, params)

	if args.o:
		for structure in structures:
			structure.write("{}_{}_{}_structure.tsv".format(args.o, structure.chrom.name.strip("chr"), tools.get_res_string(structure.chrom.res)))
	else:
		for structure in structures:
			structure.write("{}_{}_{}_structure.tsv".format(args.prefix, structure.chrom.name.strip("chr"), tools.get_res_string(structure.chrom.res)))

if __name__ == "__main__":
	main()
