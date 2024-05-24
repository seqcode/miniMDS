import numpy as np
from matplotlib import pyplot as plt
import misc

#labels = ("Chromosome3D", "mMDS", "cMDS", "miniMDS", "MOGEN", "ChromSDE")
labels = ("mMDS", "cMDS", "miniMDS", "MOGEN")
x_pos = np.arange(len(labels))

#with open("chromosome3d_chr22_10kb_time.txt") as in_file:
#	chromosomethreed_time = float(in_file.readline().strip())/60	#time in minutes
#in_file.close()

with open("mmds_chr22_10kb_time.txt") as in_file:
	mmds_time = float(in_file.readline().strip())/60	#time in minutes
in_file.close()

with open("cmds_chr22_10kb_time.txt") as in_file:
	cmds_time = float(in_file.readline().strip())/60	#time in minutes
in_file.close()

with open("minimds_chr22_10kb_time.txt") as in_file:
	minimds_time = float(in_file.readline().strip())/60	#time in minutes
in_file.close()

with open("mogen_chr22_10kb_time.txt") as in_file:
	mogen_time = float(in_file.readline().strip())/60	#time in minutes
in_file.close()

#with open("chromsde_chr22_10kb_time.txt") as in_file:
#	chromsde_time = float(in_file.readline().strip())/60	#time in minutes
#in_file.close()

#times = [chromosomethreed_time, mmds_time, cmds_time, minimds_time, mogen_time, chromsde_time]
times = [mmds_time, cmds_time, minimds_time, mogen_time]
 
#colors = ["y", "r", "g", "b", "m", "blueviolet"]
colors = ["r", "g", "b", "m"]

rects = plt.bar(x_pos, times, align="center", color = colors)
plt.yscale("log", subsy=[])
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
#plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5]), labels, fontsize=9, loc=0)
plt.legend((rects[0], rects[1], rects[2], rects[3]), labels, fontsize=9, loc=0)
plt.ylabel("Computational time (minutes)")
 
plt.savefig("Fig4.png")
