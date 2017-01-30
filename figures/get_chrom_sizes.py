import sys
sys.path.append("minimds")
import data_tools as dt

res_kb = int(sys.argv[1])
chroms = ["X", 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
with open("chrom_sizes_{}kb.txt".format(res_kb), "w") as out:
	for chrom in chroms:
		cluster = dt.clusterFromFile("data/GM12878_combined_chr{}_{}kb_cluster.tsv".format(chrom, res_kb))
		out.write(str(len(cluster.getPoints())) + "\n")
out.close()
