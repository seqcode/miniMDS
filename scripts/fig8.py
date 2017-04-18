from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
from scipy import stats as st
import misc

chroms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X")
n = len(chroms)

mmds_rs = np.zeros(n)
cmds_rs = np.zeros(n)
minimds_rs = np.zeros(n)
mogen_rs = np.zeros(n)

for i, chrom in enumerate(chroms):
	#"true" distance matrix
	bedpath = "hic_data/GM12878_combined_{}_10kb.bed".format(chrom)

	mmds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_mmds_coords.tsv".format(chrom))
	contactMat = dt.matFromBed(bedpath, mmds_cluster)
	mmds_true_mat = at.contactToDist(contactMat)
	at.makeSymmetric(mmds_true_mat)
	for j in range(len(mmds_true_mat)):	#remove diagonal
		mmds_true_mat[j,j] = 0
	mmds_distMat = misc.distMat(mmds_cluster)
	mmds_rs[i] = misc.pearson(mmds_true_mat, mmds_distMat)
	
	cmds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_cmds_coords.tsv".format(chrom))
	contactMat = dt.matFromBed(bedpath, cmds_cluster)
	cmds_true_mat = at.contactToDist(contactMat)
	at.makeSymmetric(cmds_true_mat)
	for j in range(len(cmds_true_mat)):	#remove diagonal
		cmds_true_mat[j,j] = 0
	cmds_distMat = misc.distMat(cmds_cluster)
	cmds_rs[i] = misc.pearson(cmds_true_mat, cmds_distMat)

	minimds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_minimds_coords.tsv".format(chrom))
	contactMat = dt.matFromBed(bedpath, minimds_cluster)
	minimds_true_mat = at.contactToDist(contactMat)
	at.makeSymmetric(minimds_true_mat)
	for j in range(len(minimds_true_mat)):	#remove diagonal
		minimds_true_mat[j,j] = 0
	minimds_distMat = misc.distMat(minimds_cluster)
	minimds_rs[i] = misc.pearson(minimds_true_mat, minimds_distMat)

	mogen_coords = np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_{}_10kb_rep1_coords.tsv".format(chrom))
	mogen_distMat = misc.distsFromCoords(mogen_coords)
	mogen_rs[i] = misc.pearson(mmds_true_mat, mogen_distMat)	#mMDS and MOGEN use the same matrix input procedure

chrom_sizes = np.loadtxt("chrom_sizes_10kb.txt")

fig = plt.figure()
ax = fig.add_subplot(111, frameon=False)
ax.plot(chrom_sizes, mmds_rs, linestyle="None", marker="o", markerfacecolor="r", mec="r", markersize=10, label="Standard metric MDS")
ax.plot(chrom_sizes, cmds_rs, linestyle="None", marker="o", markerfacecolor="g", mec="g", markersize=10, label="Classical MDS")
ax.plot(chrom_sizes, minimds_rs, linestyle="None", marker="o",markerfacecolor="b", mec="b", markersize=10, label="miniMDS")
ax.plot(chrom_sizes, mogen_rs, linestyle="None", marker="o",markerfacecolor="m", mec="m", markersize=10, label="MOGEN")
x_offset = 1000		#small number to prevent things from getting cut off
y_offset = 0.01
xmin = min(chrom_sizes) - x_offset
xmax = max(chrom_sizes) + x_offset
ymin = 0 - y_offset
ymax = 0.8
plt.axis([xmin, xmax, ymin, ymax])
plt.axvline(x=xmin, ymin=0, ymax=1, color="k", lw=4)
plt.axhline(y=ymin, xmin=0, xmax=1, color="k", lw=4)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=10, labelsize=14)
plt.xlabel("Number of genomic loci", fontsize=16)
plt.ylabel("Correlation between input distances and output distances", fontsize=12)
plt.legend(loc=0, numpoints=1)
plt.tight_layout()
plt.savefig("Fig8.png")
