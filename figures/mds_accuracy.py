from matplotlib import pyplot as plt
import numpy as np
import multiprocessing as mp
from sklearn import manifold
import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
import minimds as mm
import stats_tools as st

def distsFromCoords(coords):
	"""Creates distance matrix from 3D coords"""
	n = len(coords)
	distMat = np.zeros((n,n))
	for i in range(n):
		for j in range(i):
			distMat[i,j] = la.calcDistance(coords[i], coords[j])
	at.makeSymmetric(distMat)
	return distMat

def mds(distMat, metric):
	"""Performs MDS on distance matrix"""
	distMat = removeInfinite(distMat)
	mds = manifold.MDS(n_components=3, metric=metric, random_state=np.random.RandomState(seed=3), verbose=0, dissimilarity="precomputed", n_jobs=-1)
	return mds.fit(distMat).embedding_

res = int(sys.argv[1])
res_kb = res/1000
domainSmoothingParameter = 0.05
minSizeFraction = 0.01
r = 32000000
num_threads = mp.cpu_count()
args = [domainSmoothingParameter, minSizeFraction, r, num_threads]

chrom = 22

highpath = "data/GM12878_combined_{}_{}kb.bed".format(chrom, res_kb)
lowpath = "data/GM12878_combined_{}_{}kb.bed".format(chrom, res_kb*10)

cluster = dt.clusterFromBed(highpath, None, None)
contactMat = dt.matFromBed(highpath, cluster, False)
distMat = at.contactToDist(contactMat)
at.makeSymmetric(distMat)
for j in range(len(distMat)):	#remove diagonal
	distMat[j,j] = 0

cmds_distMat = distsFromCoords(st.cmds(distMat))
cmds_r = st.pearson(cmds_distMat, distMat)
mmds_distMat = distsFromCoords(mds(distMat, True))
mmds_r = st.pearson(mmds_distMat, distMat)

minimds_cluster = mm.partitionedMDS(highpath, lowpath, args)
contactMat = dt.matFromBed(highpath, minimds_cluster, False)
distMat = at.contactToDist(contactMat)
minimds_distMat = distsFromCoords(minimds_cluster.getCoords())
for j in range(len(distMat)):	#remove diagonal
	distMat[j,j] = 0
minimds_r = st.pearson(distMat, minimds_distMat)

rs = (cmds_r, mmds_r, minimds_r)

fig, ax = plt.subplots()
ax.bar(range(len(rs)), rs)
ax.set_xticklabels(("cMDS", "mMDS", "miniMDS"))
ax.set_ylabel("Pearson correlation")
plt.savefig("accuracy_{}kb.png".format(res_kb))
