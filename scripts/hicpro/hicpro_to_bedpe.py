import sys

bed_path = sys.argv[1]
mat_path = sys.argv[2]
bedpe_path = sys.argv[3]

with open(bed_path) as in_file:
    name_to_interval = {}
    for line in in_file:
        line = line.strip().split("\t")
        name_to_interval[line[3]] = "\t".join((line[0], line[1], line[2]))

with open(mat_path) as in_file, open(bedpe_path, "w") as out_file:
    for line in in_file:
        name1, name2, count = line.split("\t")
        interval1 = name_to_interval[name1]
        interval2 = name_to_interval[name2]
        out_file.write("\t".join((interval1, interval2, count)))