import sys
sys.path.append("..")
import data_tools as dt
import plotting as plot
import numpy as np

plot.plot_cluster_interactive(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_mmds_coords.tsv"), out_path="Fig9A.png")
plot.plot_cluster_interactive(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_cmds_coords.tsv"), out_path="Fig9B.png")
plot.plot_cluster_interactive(dt.clusterFromFile("hic_data/GM12878_combined_22_10kb_minimds_coords.tsv"), out_path="Fig9C.png")
plot.plot_coords_interactive(np.loadtxt("MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep1_coords.tsv"), out_path="Fig9D.png")
