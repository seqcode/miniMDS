import numpy as np
import matplotlib.pyplot as plt

def parse_time(time_string):
	split = time_string.split("m")
	mins = int(split[0])
	secs = float(split[1].split("s")[0])
	return mins + secs/60
 
labels = ('Chromosome3D', 'mMDS', 'cMDS', 'miniMDS', 'MOGEN', 'HSA', 'ChromSDE')
x_pos = np.arange(len(labels))

times = [] 
with open("chr22_100kb_times.txt".format(res_kb)) as in_file:
	for line in in_file:
		times.append(parse_time(line.strip()))
in_file.close()
 
colors = ['y', 'r', 'g', 'b', 'm', 'c', 'blueviolet']

rects = plt.bar(x_pos, times, align='center', color = colors)
plt.yscale('log', subsy=[])
plt.tick_params(top=False,bottom=False,right=False,left=False, labelbottom=False)
plt.legend((rects[0], rects[1], rects[2], rects[3], rects[4], rects[5]), labels, fontsize=9, loc=0)
plt.ylabel('Computational time (minutes)')
 
plt.savefig("Sup1.png")
