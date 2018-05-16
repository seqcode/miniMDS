from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append("..")
import data_tools as dt
import array_tools as at
import misc

bedpath = "hic_data/GM12878_combined_22_100kb.bed"

mmds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_100kb_mmds_coords.tsv")
contactMat = dt.matFromBed(bedpath, mmds_structure)
mmds_true_mat = at.contactToDist(contactMat)
at.makeSymmetric(mmds_true_mat)
for j in range(len(mmds_true_mat)):	#remove diagonal
	mmds_true_mat[j,j] = 0
mmds_distMat = misc.distMat(mmds_structure)
mmds_r = misc.pearson(mmds_true_mat, mmds_distMat)

cmds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_100kb_cmds_coords.tsv")
contactMat = dt.matFromBed(bedpath, cmds_structure)
cmds_true_mat = at.contactToDist(contactMat)
at.makeSymmetric(cmds_true_mat)
for j in range(len(cmds_true_mat)):	#remove diagonal
	cmds_true_mat[j,j] = 0
cmds_distMat = misc.distMat(cmds_structure)
cmds_r = misc.pearson(cmds_true_mat, cmds_distMat)

minimds_structure = dt.structureFromFile("hic_data/GM12878_combined_22_100kb_minimds_coords.tsv")
contactMat = dt.matFromBed(bedpath, minimds_structure)
minimds_true_mat = at.contactToDist(contactMat)
at.makeSymmetric(minimds_true_mat)
for j in range(len(minimds_true_mat)):	#remove diagonal
	minimds_true_mat[j,j] = 0
minimds_distMat = misc.distMat(minimds_structure)
minimds_r = misc.pearson(minimds_true_mat, minimds_distMat)

mogen_coords = np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_100kb_rep1_coords.tsv")
mogen_distMat = misc.distsFromCoords(mogen_coords)
mogen_r = misc.pearson(mmds_true_mat, mogen_distMat)	#mMDS and MOGEN use the same matrix input procedure

hsa_coords = np.loadtxt("hsa/GM12878_combined_22_100kb_coords.txt")
hsa_distMat = misc.distsFromCoords(hsa_coords)
hsa_r = misc.pearson(mmds_true_mat, hsa_distMat)

#chromthreed_coords = np.loadtxt("Chromosome3D/output_models/chr22_100kb/chr22_100kb_coords.tsv")
#chromthreed_distMat = misc.distsFromCoords(chromthreed_coords)
#chromthreed_r = misc.pearson(mmds_true_mat, chromthreed_distMat)

#chromsde_coords = np.loadtxt("ChromSDE/GM12878_combined_22_100kb_coords.tsv")
#chromsde_distMat = misc.distsFromCoords(chromsde_coords)
#chromsde_r = misc.pearson(mmds_true_mat, chromsde_distMat)

#labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN", "HSA", "ChromSDE")
labels = ("mMDS", "cMDS", "miniMDS", "MOGEN", "HSA")
x_pos = np.arange(len(labels))
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
