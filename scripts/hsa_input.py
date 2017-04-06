import numpy as np
import sys
sys.path.append("..")
import data_tools as dt

inpath = sys.argv[1]
outpath = sys.argv[2]

cluster = dt.clusterFromBed(inpath, None, None)
contactMat = dt.matFromBed(inpath, cluster)
maxNumDigits = int(np.ceil(np.log10(np.amax(contactMat))))
formatstring = "%" + str(maxNumDigits) + "d"
np.savetxt(outpath, contactMat, formatstring, delimiter="\t")
