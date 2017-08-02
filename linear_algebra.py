import numpy as np
import stats_tools as st

def getTransformation(cluster1, cluster2):
	"""Recovers transformation needed to align cluster1 with cluster2. Modified from http://nghiaho.com/?page_id=671"""
	pointNums1 = cluster1.getPointNums()
	pointNums2 = cluster2.getPointNums()

	intersection = [num for num in pointNums1 if num in pointNums2]

	a = []	#will hold 3D coords
	b = []	
	for num in intersection:
		a.append(cluster1.points[num-cluster1.offset].pos)
		b.append(cluster2.points[num-cluster2.offset].pos)

	a = np.mat(a)
	b = np.mat(b)

	n = a.shape[0]	#number of points

	centroid_a = np.mean(a, axis=0)
	centroid_b = np.mean(b, axis=0)

	#center the points
	aa = a - np.tile(centroid_a, (n, 1))
	bb = b - np.tile(centroid_b, (n, 1))

	h = np.transpose(aa) * bb

	u, s, vt = np.linalg.svd(h)

	r = vt.T * u.T

	# special reflection case
	#if np.linalg.det(r) < 0:
	#	vt[2,:] *= -1
	#	r = vt.T * u.T

	t = -r*centroid_a.T + centroid_b.T

	return r, t

def calcDistance(coord1, coord2):
	"""Euclidean distance between coordinates"""
	return ((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 +  (coord1[2] - coord2[2])**2)**(1./2)

def radius_of_gyration(cluster):
	coords = np.array(cluster.getCoords())
	centroid = np.mean(coords, axis=0)
	dist_sum = sum([calcDistance(coord, centroid) for coord in coords])
	return dist_sum/len(coords)
