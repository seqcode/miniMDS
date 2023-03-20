import sys

bed_path = sys.argv[1]
mat_path = sys.argv[2]
bedpe_path = sys.argv[3]

in_file = open(bed_path)
name_to_interval = {}
for line in in_file:
    chrom, chromStart, chromEnd, name = line.strip().split("\t")
    name_to_interval[name] = "\t".join((chrom, chromStart, chromEnd))
in_file.close()

in_file = open(mat_path)
out_file = open(bedpe_path, "w")
for line in in_file:
    name1, name2, count = line.split("\t")
    interval1 = name_to_interval[name1]
    interval2 = name_to_interval[name2]
    out_file.write("\t".join((interval1, interval2, count)))
in_file.close()
out_file.close()