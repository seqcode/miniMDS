***********************************************
Author: Zhang Zhizhuo (zzz2010@gmail.com)

The program is tested in Matlab2012, 
and may not support the older version of Matlab.
Comments and bug-reports are higly appreciated.
***********************************************

MainScript: ChromSDE.m

ChromSDE(trainBin,trainFreq,method_type)

Input data: (assume n loci)
trainBin: nx4 matrix, the description for each 3d point (id, chromosome, start, end)
trainFreq: nxn sparse matrix,  the normalized contact frequency matrix
method_type: 1 for quadratic SDP and 0 for linear SDP

The helper function "readpipeline_output" can generate the input data from the output(cbins, n_contact) of the Hi-C pipeline of Amos Tanay's Group
http://compgenomics.weizmann.ac.il/tanay/?page_id=283
Sample files are provided in the "data" folder.

Real data example:
[binAnno,FreqMat,normFreqMat]=readpipeline_output('data/mESC_Hind3_1000000');
ChromSDE(binAnno,normFreqMat,0)

Output data:
XXX.stat.txt: the general statistics of the predicted structures
XXX.pos:  first column is genomice locatoin = chr*10^10+position, column 2-4 are the 3D coordinates
XXX.pdb: the pdb file format for the predict 3D structure
XXX.fig: matlab figure file
XXX.png: figure for 3D structure

Note that: number of points in the pos file may differ from cbin files, as some regions are disconnected to other regions, and ChromSDE only consider the largest connected component defined by the contact matrix.


Simuation Demo:
SimulationDemo.m



