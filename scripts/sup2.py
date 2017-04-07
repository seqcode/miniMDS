import numpy as np
from matplotlib import pyplot as plt

with open("chromosome3d_chr22_100kb_memory.txt") as in_file:
	chromthreed_mem = float(in_file.readline().strip())
in_file.close()

with open("mmds_chr22_100kb_memory.txt") as in_file:
	mmds_mem = float(in_file.readline().strip())
in_file.close()

with open("cmds_chr22_100kb_memory.txt") as in_file:
	cmds_mem = float(in_file.readline().strip())
in_file.close()

with open("minimds_chr22_100kb_memory.txt") as in_file:
	minimds_mem = float(in_file.readline().strip())
in_file.close()

with open("mogen_chr22_100kb_memory.txt") as in_file:
	mogen_mem = float(in_file.readline().strip())
in_file.close()

#with open("chromsde_chr22_100kb_memory.txt") as in_file:
#	chromsde_mem = float(in_file.readline().strip())
#in_file.close()
 
#labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN", "ChromSDE")
labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN")
x_pos = np.arange(len(objects))
#memory = [chromthreed_mem, mmds_mem, cmds_mem, minimds_mem, mogen_mem, chromsde_mem] 
memory = [chromthreed_mem, mmds_mem, cmds_mem, minimds_mem, mogen_mem]
 
#colors = ["y", "r", "g", "b", "m", "blueviolet"]
colors = ["y", "r", "g", "b", "m"]

rects = plt.bar(x_pos, memory, align="center", color = colors)
plt.title("100-Kbp resolution", fontsize=12)
plt.yscale("log", subsy=[])
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
#plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5]), labels, fontsize=9, loc=0)
plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4]), labels, fontsize=9, loc=0)
plt.ylabel("Memory usage (Mb)")
plt.savefig("Sup2.png")
