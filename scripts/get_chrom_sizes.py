import sys
sys.path.append("..")
import data_tools as dt

res_kb = int(sys.argv[1])
with open("chrom_sizes_{}kb.txt".format(res_kb), "w") as out:
	for chrom in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]:
		cluster = dt.clusterFromFile("hic_data/GM12878_combined_chr{}_{}kb_cluster.tsv".format(chrom, res_kb))
		out.write(str(len(cluster.getPoints())) + "\n")
out.close()
