import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot
import numpy as np
import misc

mmds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_mmds_coords.tsv")
cmds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_cmds_coords.tsv")
minimds_cluster = dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_minimds_coords.tsv")

mmds_res = mmds_cluster.chrom.res
cmds_res = cmds_cluster.chrom.res
minimds_res = minimds_cluster.chrom.res

assert mmds_res == cmds_res == minimds_res

res = mmds_res

plot.plot_cluster_interactive(mmds_cluster, out_path="Fig9A.png")
plot.plot_cluster_interactive(cmds_cluster, out_path="Fig9B.png")
plot.plot_cluster_interactive(minimds_cluster, out_path="Fig9C.png")
misc.plot_coords_interactive(np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep1_coords.tsv"), res, out_path="Fig9D.png")
