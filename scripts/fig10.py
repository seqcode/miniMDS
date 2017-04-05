import sys
sys.path.append("..")
import plotting as plot
import data_tools as dt

chroms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
clusters = [dt.clusterFromFile("hic_data/GM12878_combined_{}_10kb_cluster.tsv".format(chrom)) for chrom in chroms]
plot.plot_clusters_interactive(clusters)
