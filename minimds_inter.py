import numpy as np
import data_tools as dt
import linear_algebra as la
import sys
from sklearn import manifold
import tools
import argparse
import minimds as mm

def get_inter_mat(intra_prefix, inter_prefix, res, clusters, offsets):
	res_string = tools.get_res_string(res)

	names = [cluster.chrom.name for cluster in clusters]
	n = len(names)
	for i in range(n):
		if names[i].startswith("chr"):
			names[i] = names[i][3:len(names[i])]	#remove "chr"

	#fill matrix
	total_len = sum([len(cluster.getPoints()) for cluster in clusters])
	mat = np.zeros((total_len, total_len))
	for i in range(n):
		for j in range(i+1):
			if i == j:
				path = "{}_{}_{}.bed".format(intra_prefix, names[i], res_string)
			else:
				path = "{}_{}_{}_{}.bed".format(inter_prefix, names[j], names[i], res_string)
			print "Reading {}".format(path)
			with open(path) as bed:
				for line in bed:
					line = line.strip().split()
					loc1 = int(line[4])
					loc2 = int(line[1])
					index1 = clusters[i].getIndex(loc1)
					index2 = clusters[j].getIndex(loc2)
					if index1 is not None and index2 is not None:
						index1 += offsets[i]
						index2 += offsets[j]
						row = index1
						col = index2
						mat[row, col] += float(line[6])
			bed.close()
	return mat

def interMDS(names, inter_prefix, intra_prefix, inter_res, intra_res, intra_low_res=None, args=None):
	inter_res_string = tools.get_res_string(inter_res)
	intra_res_string = tools.get_res_string(intra_res)
	if intra_low_res is None:
		intra_low_res_string = None
	else:
		intra_low_res_string = tools.get_res_string(intra_low_res)

	#get low-res clusters from intra files
	low_clusters = [dt.clusterFromBed("{}_{}_{}.bed".format(intra_prefix, name, inter_res_string), None, None) for name in names]

	#for correct indexing
	n = len(names)
	offsets = np.zeros(n, dtype=np.int32)
	for i in range(1, n):
		offsets[i] = offsets[i-1] + len(low_clusters[i-1].getPoints())

	inter_mat = get_inter_mat(intra_prefix, inter_prefix, inter_res, low_clusters, offsets)

	#perform MDS at low resolution on all chroms
	mm.infer_clusters(inter_mat, low_clusters, offsets)

	#perform MDS at high resolution on each chrom
	high_clusters = []
	inferred_low_clusters = []
	ts = []
	for true_low, name in zip(low_clusters, names):
		path = "{}_{}_{}.bed".format(intra_prefix, name, intra_res_string)
		if intra_low_res_string is None:
			high_cluster = mm.fullMDS(path, False)
		else:
			low_path = "{}_{}_{}.bed".format(intra_prefix, name, intra_low_res_string)
			high_cluster = mm.partitionedMDS(path, low_path, args)
		high_clusters.append(high_cluster)
		inferred_low = dt.highToLow(high_cluster, true_low.chrom.res/high_cluster.chrom.res)
		inferred_low_clusters.append(inferred_low)

		#rescale
		rescaling_factor = la.radius_of_gyration(true_low)/la.radius_of_gyration(inferred_low)
		rescaled_coords = [rescaling_factor * coord for coord in inferred_low.getCoords()]
		for i, point in enumerate(inferred_low.getPoints()):
			point.pos = rescaled_coords[i]

		r, t, reflect = la.getTransformation(inferred_low, true_low)
		high_cluster.transform(r, None, reflect)	#do not translate now (need to rescale)
		ts.append(t)	

	#translate (with rescaling)
	low_rgs = np.array([la.radius_of_gyration(cluster) for cluster in low_clusters])
	high_rgs = np.array([la.radius_of_gyration(cluster) for cluster in high_clusters])
	scaling_factor = np.mean(high_rgs/low_rgs)
	for high_cluster, t in zip(high_clusters, ts):
		high_cluster.transform(None, scaling_factor*t, False)	#rescale translation

	return high_clusters

def main():
	parser = argparse.ArgumentParser(description="Reconstruct 3D coordinates from normalized interchromosomal Hi-C BED files.")
	parser.add_argument("inter_prefix", help="prefix of interchromosomal Hi-C BED file")
	parser.add_argument("intra_prefix", help="prefix of intrachromosomal Hi-C BED file")
	parser.add_argument("inter_res", type=int, help="resolution of interchromosomal BED files (bp)")
	parser.add_argument("intra_res", type=int, help="resolution of intrachromosomal BED files (bp)")
	parser.add_argument("-c", action="append", default=[], help="Names of chromosomes to use, e.g. 1 (default: all human chromosomes other than Y)")
	parser.add_argument("-l", type=int, help="resolution of low-res intrachromosomal files (bp) (for partitioned MDS only)")
	parser.add_argument("-p", type=float, default=0.1, help="domain size parameter: larger value means fewer clusters created (for partitioned MDS only)")
	parser.add_argument("-m", type=float, default=0.05, help="minimum domain size parameter: prevents clusters from being too small (for partitioned MDS only)")
	parser.add_argument("-o", help="prefix of output file")
	parser.add_argument("-r", default=32000000, help="maximum RAM to use (in kb)")
	parser.add_argument("-n", default=3, help="Number of threads")
	args = parser.parse_args()

	if len(args.c) == 0:
		chrom_names = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
	else:
		chrom_names = args.c

	if args.l is None:	#not partitioned
		clusters = interMDS(chrom_names, args.inter_prefix, args.intra_prefix, args.inter_res, args.intra_res)
	else:	#partitioned
		params = (args.p, args.m, args.r, args.n)
		names = ("Domain size parameter", "Minimum domain size", "Maximum memory", "Number of threads")
		intervals = ((0,1), (0,1), (0, None), (0, None))
		if not tools.args_are_valid(params, names, intervals):
			sys.exit(0)

		clusters = interMDS(chrom_names, args.inter_prefix, args.intra_prefix, args.inter_res, args.intra_res, args.l, params)

	if args.o is not None:
		for cluster in clusters:
			cluster.write("{}_{}_{}_cluster.tsv".format(args.o, cluster.chrom.name, tools.get_res_string(cluster.chrom.res)))

if __name__ == "__main__":
	main()
