import numpy as np

def getTransformation(structure1, structure2):
	"""Recovers transformation needed to align structure1 with structure2. Modified from http://nghiaho.com/?page_id=671"""
	pointNums1 = structure1.getPointNums()
	pointNums2 = structure2.getPointNums()

	intersection = [num for num in pointNums1 if num in pointNums2]

	a = []	#will hold 3D coords
	b = []	
	for num in intersection:
		a.append(structure1.points[num-structure1.offset].pos)
		b.append(structure2.points[num-structure2.offset].pos)

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

	t = -r*centroid_a.T + centroid_b.T

	return r, t

def calcDistance(coord1, coord2):
	"""Euclidean distance between coordinates"""
	return ((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 +  (coord1[2] - coord2[2])**2)**(1./2)

def radius_of_gyration(structure):
	coords = np.array(structure.getCoords())
	centroid = np.mean(coords, axis=0)
	dist_sum = sum([calcDistance(coord, centroid) for coord in coords])
	return dist_sum/len(coords)

def cmds(distMat):
	"""Modified from http://www.nervouscomputer.com/hfs/cmdscale-in-python/"""
	# Number of points                                                                        
	n = len(distMat)

	# Centering matrix                                                                        
	h = np.eye(n) - np.ones((n, n))/n

	# YY^T                                                                                    
	b = -h.dot(distMat**2).dot(h)/2

	# Diagonalize                                                                             
	evals, evecs = np.linalg.eigh(b)

	# Sort by eigenvalue in descending order                                                  
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	evecs = evecs[:,idx]

	return np.array([evecs[:,0]*evals[0]**(1./2), evecs[:,1]*evals[1]**(1./2), evecs[:,2]*evals[2]**(1./2)]).T

def change_coordinate_system(n, coords):
	"""Rotate 3-D coords such that vector n aligns to z-axis""" 
	origin = [0,0,0]
	len_n = calcDistance(n, origin)
	n_hat = n/len_n
	i_hat = [1,0,0]		#x-axis
	j_hat = [0,1,0]		#y-axis
	k_hat = [0,0,1]		#z-axis
	theta = np.arccos(np.dot(k_hat, n_hat))
	b = np.cross(k_hat, n_hat)
	len_b = calcDistance(b, origin)
	b_hat = b/len_b

	q0 = np.cos(theta/2)
	q1 = np.sin(theta/2)*b_hat[0]
	q2 = np.sin(theta/2)*b_hat[1]
	q3 = np.sin(theta/2)*b_hat[2]

	Q = np.matrix([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)], [2*(q2*q1 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 - q0*q1)], [2*(q3*q1 - q0*q2), 2*(q3*q2 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])

	i_hat = np.matrix(i_hat).transpose()
	j_hat = np.matrix(j_hat).transpose()
	k_hat = np.matrix(k_hat).transpose()

	u_hat = Q*i_hat
	v_hat = Q*j_hat
	w_hat = Q*k_hat

	new_coords = np.zeros_like(coords)
	for i in range(len(coords)):
		p = coords[i]
		new_coords[i] = [np.dot(p, u_hat), np.dot(p, v_hat), np.dot(p, w_hat)]
	return new_coords
