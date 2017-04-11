import sys
sys.path.append("..")
import data_tools as dt
import misc

in_path = sys.argv[1]
out_path = sys.argv[2]

chrom = dt.intraChromFromBed(in_path, None)

with open(in_path) as in_file:
	with open(out_path, "w") as out_file:
		for line in in_file:
			line = line.strip().split()
			loc1 = int(line[1])
			loc2 = int(line[4])
			num1 = misc.getPointNum(chrom, loc1)
			num2 = misc.getPointNum(chrom, loc2)
			out_file.write("\t".join((str(num1), str(num2), line[6])) + "\n")
	out_file.close()
in_file.close()
