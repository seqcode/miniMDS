import sys
import numpy as np
from sklearn import manifold
import argparse
import data_tools as dt
import array_tools as at
import tad
import linear_algebra as la
import tools
import stats_tools as st
import pymp
import multiprocessing as mp

def infer_clusters(contactMat, clusters, offsets, alpha, classical=False):
	"""Infers 3D coordinates for multiple clusters with same contact matrix"""
	assert sum([len(cluster.getPointNums()) for cluster in clusters]) == len(contactMat)

	at.makeSymmetric(contactMat)
	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)
	at.makeSymmetric(distMat)

	if classical:	#classical MDS
		coords = st.cmds(distMat)
	else:
		mds = manifold.MDS(n_components=3, metric=True, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=-1)
		coords = mds.fit_transform(distMat)

	for offset, cluster in zip(offsets, clusters):
		for i in range(len(cluster.getPoints())):	
			cluster.getPoints()[i].pos = coords[i + offset]

def infer_cluster(contactMat, cluster, alpha, classical=False):
	"""Infers 3D coordinates for one cluster"""
	assert len(cluster.getPointNums()) == len(contactMat)

	at.makeSymmetric(contactMat)
	rowsums = np.array([sum(row) for row in contactMat])
	assert len(np.where(rowsums == 0)[0]) == 0 

	distMat = at.contactToDist(contactMat, alpha)
	at.makeSymmetric(distMat)

	if classical:	#classical MDS
		coords = st.cmds(distMat)
	else:
		mds = manifold.MDS(n_components=3, metric=True, random_state=np.random.RandomState(), verbose=0, dissimilarity="precomputed", n_jobs=-1)
		coords = mds.fit_transform(distMat)

	for i in range(len(cluster.getPoints())):	
		cluster.getPoints()[i].pos = coords[i]

def fullMDS(path, classical, alpha):
	"""MDS without partitioning"""
	cluster = dt.clusterFromBed(path, None, None)
	contactMat = dt.matFromBed(path, cluster)
	infer_cluster(contactMat, cluster, alpha, classical)
	return cluster
	
def partitionedMDS(path, args):
	"""Partitions cluster into subclusters and performs MDS"""
	domainSmoothingParameter = args[0]
	minSizeFraction = args[1]
	maxmemory = args[2]
	num_threads = args[3]
	alpha = args[4]
	res_ratio = args[5]

	#create low-res cluster
	low_chrom = dt.chromFromBed(path)
	low_chrom.res *= res_ratio
	lowCluster = dt.clusterFromBed(path, low_chrom, None)

	#get TADs
	low_contactMat = dt.matFromBed(lowpath, lowCluster)
	lowTads = tad.getDomains(low_contactMat, lowCluster, domainSmoothingParameter, minSizeFraction)		#low subclusters

	#create high-res chrom
	size, res = dt.basicParamsFromBed(path)
	highChrom = dt.ChromParameters(lowCluster.chrom.minPos, lowCluster.chrom.maxPos, res, lowCluster.chrom.name, size)

	#create high-res cluster
	res_ratio = lowCluster.chrom.res/highChrom.res
	highTads = lowTads * res_ratio
	highCluster = dt.clusterFromBed(path, highChrom, highTads)

	#create compatible subclusters
	tad.subclustersFromTads(highCluster, lowCluster, lowTads)

	infer_cluster(low_contactMat, lowCluster, alpha)
	print "Low-resolution MDS complete"

	highSubclusters = pymp.shared.list(highCluster.clusters)
	lowSubclusters = pymp.shared.list(lowCluster.clusters)

	numSubclusters = len(highCluster.clusters)
	num_threads = min((num_threads, mp.cpu_count(), numSubclusters))	#don't exceed number of requested threads, available threads, or clusters
	with pymp.Parallel(num_threads) as p:
		for subclusternum in p.range(numSubclusters):
			highSubcluster = highSubclusters[subclusternum]
			trueLow = lowSubclusters[subclusternum]

			#perform MDS individually
			cluster_contactMat = dt.matFromBed(path, highSubcluster)	#contact matrix for this cluster only
			infer_cluster(cluster_contactMat, highSubcluster, 2)

			#approximate as low resolution
			inferredLow = dt.highToLow(highSubcluster, res_ratio)

			#rescale
			scaling_factor = la.radius_of_gyration(trueLow)/la.radius_of_gyration(inferredLow)
			for i, point in enumerate(inferredLow.points):
				if point != 0:
					x, y, z = point.pos
					inferredLow.points[i].pos = (x*scaling_factor, y*scaling_factor, z*scaling_factor)
	
			#recover the transformation for inferred from true low cluster
			r, t = la.getTransformation(inferredLow, trueLow)
			t *= res_ratio**(2./3)	#rescale

			#transform high cluster
			highSubcluster.transform(r, t)
			highSubclusters[subclusternum] = highSubcluster

			print "MDS performed on cluster {} of {}".format(subclusternum + 1, numSubclusters)

	highCluster.setClusters(highSubclusters)

	return highCluster

def main():
	parser = argparse.ArgumentParser(description="Reconstruct 3D coordinates from normalized intrachromosomal Hi-C BED files.")
	parser.add_argument("path", help="path to intrachromosomal Hi-C BED file")
	parser.add_argument("--classical", action="store_true", help="use classical MDS (default: metric MDS)")
	parser.add_argument("--full", action="store_true", help="use full MDS (default: partitioned MDS)")
	parser.add_argument("-l", help="low resolution/high resolution")
	parser.add_argument("-p", type=float, default=0.1, help="domain size parameter: larger value means fewer clusters created (for partitioned MDS only)")
	parser.add_argument("-m", type=float, default=0.05, help="minimum domain size parameter: prevents clusters from being too small (for partitioned MDS only)")
	parser.add_argument("-o", help="path to output file")
	parser.add_argument("-r", default=32000000, help="maximum RAM to use (in kb)")
	parser.add_argument("-n", default=3, help="number of threads")
	parser.add_argument("-a", type=float, default=4, help="alpha factor for converting contact frequencies to physical distances")
	args = parser.parse_args()

	if args.full:	#not partitioned
		cluster = fullMDS(args.path, args.classical, args.a)

	else:	#partitioned
		params = (args.p, args.m, args.r, args.n, args.a, args.l)
		names = ("Domain size parameter", "Minimum domain size", "Maximum memory", "Number of threads", "Alpha", "Resolution ratio")
		intervals = ((0, 1), (0, 1), (0, None), (0, None), (1, None), (1, None))
		if not tools.args_are_valid(params, names, intervals):
			sys.exit(0)

		cluster = partitionedMDS(args.path, params)
	
	if args.o is not None:
		cluster.write(args.o)

if __name__ == "__main__":
	main()
