import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

in_path = sys.argv[1]
out_path = sys.argv[2]

contactMat = dt.matFromBed(in_path)
np.savetxt(out_path, contactMat, delimiter="\t")
