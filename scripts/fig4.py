import numpy as np
from matplotlib import pyplot as plt
import misc

labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN", "ChromSDE")
x_pos = np.arange(len(labels))

times = [] 
with open("chr22_10kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		times.append(misc.parse_time(line.strip()))
in_file.close()
 
colors = ["y", "r", "g", "b", "m", "blueviolet"]

rects = plt.bar(x_pos, times, align="center", color = colors)
plt.yscale("log", subsy=[])
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5]), labels, fontsize=9, loc=0)
plt.ylabel("Computational time (minutes)")
 
plt.savefig("Fig4.png")
