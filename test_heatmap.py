import sys
sys.path.append("../figures")
import heatmap as hm
import numpy as np

mat = np.load("21_22_inter_mat.npy")
hm.heatMapFromMat(mat, 5000, None, None)
