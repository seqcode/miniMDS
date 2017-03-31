import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

in_path = sys.argv[1]
out_path = sys.argv[2]

cluster = dt.clusterFromBed(in_path, None, None)
contactMat = dt.matFromBed(in_path, cluster)
np.savetxt(out_path, contactMat, delimiter="\t")
