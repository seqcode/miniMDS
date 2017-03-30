import numpy as np
from matplotlib import pyplot as plt

res_kb = 10

chrom_sizes = np.loadtxt("chrom_sizes_{}kb.txt".format(res_kb))

mmds_memories = np.loadtxt("mmds_{}kb_memory.txt".format(res_kb))/10**6
cmds_memories = np.loadtxt("cmds_{}kb_memory.txt".format(res_kb))/10**6
minimds_memories = np.loadtxt("minimds_{}kb_memory.txt".format(res_kb))/10**6
mogen_memories = np.loadtxt("mogen_{}kb_memory.txt".format(res_kb))/10**6

fig = plt.figure()
ax = fig.add_subplot(111, frameon=False)
ax.plot(chrom_sizes, mmds_memories, linestyle="None", marker="o", markerfacecolor="r", mec="r", markersize=10, label="Standard metric MDS")
ax.plot(chrom_sizes, cmds_memories, linestyle="None", marker="o", markerfacecolor="g", mec="g", markersize=10, label="Classical MDS")
ax.plot(chrom_sizes, minimds_memories, linestyle="None", marker="o",markerfacecolor="b", mec="b", markersize=10, label="miniMDS")
ax.plot(chrom_sizes, mogen_memories, linestyle="None", marker="o",markerfacecolor="m", mec="m", markersize=10, label="MOGEN")
x_offset = 300		#small number to prevent things from getting cut off
y_offset = 3
xmin = min(chrom_sizes) - x_offset
xmax = max(chrom_sizes) + x_offset
ymin = 0 - y_offset
ymax = max((max(mmds_memories), max(cmds_memories), max(minimds_memories))) + y_offset
plt.axis([xmin, xmax, ymin, ymax])
plt.axvline(x=xmin, ymin=0, ymax=1, color="k", lw=4)
plt.axhline(y=ymin, xmin=0, xmax=1, color="k", lw=4)
plt.tick_params(direction="out", top=False, right=False, length=12,width=3, pad=10, labelsize=14)
plt.xlabel("Number of genomic loci", fontsize=16)
plt.ylabel("Computational memory (Gb)", fontsize=16)
plt.legend(loc=0, numpoints=1)
plt.tight_layout()
plt.savefig("memory_{}kb.png".format(res_kb))
