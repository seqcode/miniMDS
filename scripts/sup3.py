from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
import misc

bedpath = "hic_data/GM12878_combined_22_10kb.bed"

#"true" distance matrix
cluster = dt.clusterFromBed(bedpath, None, None)
contactMat = dt.matFromBed(bedpath, cluster)
distMat = at.contactToDist(contactMat)
at.makeSymmetric(distMat)
for j in range(len(distMat)):	#remove diagonal
	distMat[j,j] = 0

#chromthreed_coords = np.loadtxt("Chromosome3D/output_models/chr22_100kb/chr22_100kb_coords.tsv")
#chromthreed_distMat = misc.distsFromCoords(chromthreed_coords)
#chromthreed_r = misc.pearson(distMat, chromthreed_distMat)

mmds_distMat = misc.distMat(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_mmds_coords.tsv"))
mmds_r = misc.pearson(distMat, mmds_distMat)

cmds_distMat = misc.distMat(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_cmds_coords.tsv"))
cmds_r = misc.pearson(distMat, cmds_distMat)

minimds_distMat = misc.distMat(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_minimds_coords.tsv"))
minimds_r = misc.pearson(distMat, minimds_distMat)

mogen_coords = np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep1_coords.tsv")
mogen_distMat = misc.distsFromCoords(mogen_coords)
mogen_r = misc.pearson(distMat, mogen_distMat)

hsa_coords = np.loadtxt("hsa/GM12878_combined_22_100kb_coords.txt")
hsa_distMat = misc.distsFromCoords(hsa_coords)
hsa_r = misc.pearson(distMat, hsa_distMat)

#chromsde_coords = np.loadtxt("ChromSDE/GM12878_combined_22_100kb_coords.tsv")
#chromsde_distMat = misc.distsFromCoords(chromsde_coords)
#chromsde_r = misc.pearson(distMat, chromsde_distMat)

#labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN", "HSA", "ChromSDE")
labels = ("mMDS", "cMDS", "miniMDS", "MOGEN", "HSA")
x_pos = np.arange(len(objects))
#rs = [chromthreed_r, mmds_r, cmds_r, minimds_r, mogen_r, hsa_r, chromsde_r] 
rs = [mmds_r, cmds_r, minimds_r, mogen_r, hsa_r] 
#colors = ["y", "r", "g", "b", "c", "m", "blueviolet"]
colors = ["r", "g", "b", "c", "m"]
rects = plt.bar(x_pos, rs, align="center", color = colors)
plt.title("100-Kbp resolution", fontsize=12)
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
#plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5], rects[6]), labels, fontsize=9, loc=0)
plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4]), labels, fontsize=9, loc=0)
plt.ylabel("Correlation between input distances and output distances")
plt.savefig("Sup3.png")
