from matplotlib import pyplot as plt
from scipy import stats as st
import numpy as np
import sys
sys.path.append("../minimds")
import array_tools as at
import linear_algebra as la
import misc

def simCoords(n):
	"Simulates n loci in 3-D space with ideal chain model"
	coords=np.zeros((n,3))
	for i in range(n):
		change = np.random.uniform(-10, 10, 3)
		if i==0:
			prevCoords=np.asarray([0,0,0])
		else:
			prevCoords=coords[i-1]
		coords[i]=prevCoords+change
	return coords

def distToContact(distMat, alpha=-3):
	"""Convert distance matrix to contact matrix. Matrix must be symmetric."""
	contactMat = np.zeros_like(distMat)
	n1 = len(contactMat)
	n2 = len(contactMat[0])
	for i in range(n1):
		for j in range(n2):
			if distMat[i,j] != 0:
				contactMat[i,j] = distMat[i,j]**alpha		#see Varoquaux, et al (2014)
	return contactMat

def simDistMat(n):
	coords = simCoords(n)
	distMat = np.zeros((n, n))
	for i in range(n):
		for j in range(i):
			distMat[i,j] = la.calcDistance(coords[i], coords[j])
	return distMat

def pearson(mat1, mat2):
	"""Pearson correlation between elements of matrices, ignoring zeroes"""
	at.makeSymmetric(mat1)
	at.makeSymmetric(mat2)
	vector1 = np.ravel(mat1)
	vector2 = np.ravel(mat2)
	n = len(vector1)
	assert n == len(vector2)
	nonzero1 = np.where(vector1 != 0)[0]
	nonzero2 = np.where(vector2 != 0)[0]
	intersection = list(set(nonzero1) & set(nonzero2))
	vector1_trimmed = [vector1[i] for i in intersection]
	vector2_trimmed = [vector2[i] for i in intersection]
	return st.linregress(vector1_trimmed, vector2_trimmed)

alphas = [-4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5]
m = 500
n = len(alphas)
mmds_rs = np.zeros(n)
nmds_rs = np.zeros(n)
for i in range(n):
	distMat = simDistMat(m)
	contactMat = distToContact(distMat, alphas[i])
	inferred_distMat = at.contactToDist(contactMat)
	at.makeSymmetric(inferred_distMat)
	dissimMat = at.contactToDist(contactMat, -1)
	at.makeSymmetric(dissimMat)
	mmdsDists = misc.distsFromCoords(misc.mds(inferred_distMat, True))
	nmdsDists = misc.distsFromCoords(misc.mds(dissimMat, False))
	slope, intercept, mmds_r, p, se = pearson(mmdsDists, distMat)
	slope, intercept, nmds_r, p, se = pearson(nmdsDists, distMat)
	mmds_rs[i] = mmds_r
	nmds_rs[i] = nmds_r

fig = plt.figure()
ax = fig.add_subplot(111, frameon=False)
ax.plot(alphas, mmds_rs, linestyle="None", marker="o", markerfacecolor="r", mec="r", markersize=10, label="Standard metric MDS")
ax.plot(alphas, nmds_rs, linestyle="None", marker="o", markerfacecolor=(0.5,0,1), mec=(0.5,0,1), markersize=10, label="Classical MDS")
x_offset = 0.1	#small number to prevent things from getting cut off
y_offset = 0.05
xmin = min(alphas) - x_offset
xmax = max(alphas) + x_offset
ymin = -0.2 - y_offset
ymax = 1 + y_offset
plt.axis([xmin, xmax, ymin, ymax])
plt.axvline(x=xmin, ymin=0, ymax=1, color="k", lw=4)
plt.axhline(y=ymin, xmin=0, xmax=1, color="k", lw=4)
plt.tick_params(direction="out", top=False, right=False, length=12,width=3, pad=10, labelsize=14)
plt.xlabel("Alpha", fontsize=16)
plt.ylabel("Pearson correlation", fontsize=16)
plt.legend(loc=6, numpoints=1)
plt.tight_layout()
plt.savefig("alphas.png")
plt.show()
