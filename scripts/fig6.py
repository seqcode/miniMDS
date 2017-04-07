import numpy as np
from matplotlib import pyplot as plt

res_kb = 10

chrom_sizes = np.loadtxt("chrom_sizes_{}kb.txt".format(res_kb))

mmds_times = []
with open("mmds_{}kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		mmds_times.append(float(line.strip())/60)
in_file.close()

cmds_times = []
with open ("cmds_{}kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		cmds_times.append(float(line.strip())/60)
in_file.close()

minimds_times = []
with open("minimds_{}kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		minimds_times.append(float(line.strip())/60)
in_file.close()

mogen_times = []
with open("mogen_{}kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		mogen_times.append(float(line.strip())/60)
in_file.close()

fig = plt.figure()
ax = fig.add_subplot(111, frameon=False)
ax.plot(chrom_sizes, mmds_times, linestyle="None", marker="o", markerfacecolor="r", mec="r", markersize=10, label="Standard metric MDS")
ax.plot(chrom_sizes, cmds_times, linestyle="None", marker="o", markerfacecolor="g", mec="g", markersize=10, label="Classical MDS")
ax.plot(chrom_sizes, minimds_times, linestyle="None", marker="o", markerfacecolor="b", mec="b", markersize=10, label="miniMDS")
ax.plot(chrom_sizes, mogen_times, linestyle="None", marker="o", markerfacecolor="m", mec="m", markersize=10, label="MOGEN")
x_offset = 1000		#small number to prevent things from getting cut off
y_offset = 5
xmin = min(chrom_sizes) - x_offset
xmax = max(chrom_sizes) + x_offset
ymin = 0 - y_offset
ymax = max((max(mmds_times), max(cmds_times), max(minimds_times))) + y_offset
plt.axis([xmin, xmax, ymin, ymax])
plt.axvline(x=xmin, ymin=0, ymax=1, color="k", lw=4)
plt.axhline(y=ymin, xmin=0, xmax=1, color="k", lw=4)
plt.tick_params(direction="out", top=False, right=False, length=12, width=3, pad=10, labelsize=14)
plt.xlabel("Number of genomic loci", fontsize=16)
plt.ylabel("Time (minutes)", fontsize=16)
plt.legend(loc=0, numpoints=1)
plt.tight_layout()
plt.savefig("Fig6.png".format(res_kb))
