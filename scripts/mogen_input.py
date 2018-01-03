import sys
sys.path.append("..")
import data_tools as dt

in_path = sys.argv[1]
out_path = sys.argv[2]

chrom = dt.chromFromBed(in_path)

with open(in_path) as in_file:
	with open(out_path, "w") as out_file:
		for line in in_file:
			line = line.strip().split()
			loc1 = int(line[1])
			loc2 = int(line[4])
			num1 = chrom.getPointNum(loc1)
			num2 = chrom.getPointNum(loc2)
			out_file.write("\t".join((str(num1), str(num2), line[6])) + "\n")
	out_file.close()
in_file.close()
