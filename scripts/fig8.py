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
	cluster = dt.clusterFromBed(bedpath, None, None)
	contactMat = dt.matFromBed(bedpath, cluster)
	distMat = at.contactToDist(contactMat)
	at.makeSymmetric(distMat)
	for j in range(len(distMat)):	#remove diagonal
		distMat[j,j] = 0

	mmds_distMat = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_mmds_coords.tsv".format(chrom)).distMat()
	mmds_rs[i] = misc.pearson(distMat, mmds_distMat)
	
	cmds_distMat = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_cmds_coords.tsv".format(chrom)).distMat()
	cmds_rs[i] = misc.pearson(distMat, cmds_distMat)

	minimds_distMat = dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_minimds_coords.tsv".format(chrom)).distMat()
	minimds_rs[i] = misc.pearson(distMat, minimds_distMat)

	mogen_distMat = misc.distsFromCoords("MOGEN/examples/hiC/output/GM12878_combined_{}_10kb_rep1_coords.tsv".format(chrom))
	mogen_rs[i] = misc.pearson(distMat, mogen_distMat)

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
