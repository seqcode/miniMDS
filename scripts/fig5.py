from matplotlib import pyplot as plt
import misc

def rep_correlation(coords1, coords2):
	dists1 = misc.distsFromCoords(coords1)
	dists2 = misc.distsFromCoords(coords2)

	assert dists1.shape == dists2.shape

	#convert to vectors
	dists1 = dists1.flatten()
	dists2 = dists2.flatten()

	#remove zeroes
	indices = np.where(dists1 != 0)[0]
	dists1 = dists1[indices]
	dists2 = dists2[indices]

	assert len(np.where(dists2 == 0)[0]) == 0

	r, p = st.pearsonr(dists1, dists2)
	return r


#labels = ("Chromosome3D", "mMDS", "miniMDS", "MOGEN", "HSA", "ChromSDE")
labels = ("mMDS", "miniMDS", "MOGEN", "HSA")
n = len(labels)
rs = np.zeros(n)

#Chromosome3D
#coords1 = np.loadtxt("Chromosome3D/output_models/chr22_10kb_rep1/rep1_coords.tsv")
#coords2 = np.loadtxt("Chromosome3D/output_models/chr22_10kb_rep1/rep2_coords.tsv")
#rs[0] = rep_correlation(coords1, coords2)

#mMDS
coords1 = np.loadtxt("hic_data/GM12878_combined_22_10kb_mmds_rep1.tsv")
coords2 = np.loadtxt("hic_data/GM12878_combined_22_10kb_mmds_rep2.tsv")
rs[1] = rep_correlation(coords1, coords2)

#miniMDS
coords1 = np.loadtxt("hic_data/GM12878_combined_22_10kb_minimds_rep1.tsv")
coords2 = np.loadtxt("hic_data/GM12878_combined_22_10kb_minimds_rep2.tsv")
rs[2] = rep_correlation(coords1, coords2)

#MOGEN
coords1 = np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep1_coords.tsv")
coords2 = np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep2_coords.tsv")
rs[3] = rep_correlation(coords1, coords2)

#HSA
coords1 = np.loadtxt("hsa/GM12878_combined_22_10kb_rep1_coords.txt")
coords2 = np.loadtxt("hsa/GM12878_combined_22_10kb_rep2_coords.txt")
rs[4] = rep_correlation(coords1, coords2)

#ChromSDE
#coords1 = np.loadtxt("ChromSDE/GM12878_combined_22_10kb_rep1_coords.tsv")
#coords2 = np.loadtxt("ChromSDE/GM12878_combined_22_10kb_rep2_coords.tsv")
#rs[5] = rep_correlation(coords1, coords2)

x_pos = range(n) 
colors = ["y", "r", "b", "c", "m", "blueviolet"]
rects = plt.bar(x_pos, rs, align="center", color = colors)
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
#plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5]), labels, fontsize=8, loc=3)
plt.legend((rects[0], rects[1], rects[2], rects[3]), labels, fontsize=8, loc=3)
plt.ylabel("Correlation between iterations")
 
plt.savefig("Fig5.png")
