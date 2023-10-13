import sys
sys.path.append("..")
import numpy as np
import argparse
import tools
import os

def get_chrom_num(chrom):
	if chrom == "X":
		return 23
	else:
		return int(chrom)

def normalize(chrom1, chrom2, rawpath, krpath1, krpath2, res, outpath):
	kr1 = np.loadtxt(krpath1)
	if krpath2 is None:
		kr2 = kr1
	else:
		kr2 = np.loadtxt(krpath2)
	with open(rawpath) as raw:
		with open(outpath, "w") as out:
			for line in raw:
				line = line.split()
				loc1 = line[0]
				loc2 = line[1]
				norm1 = kr1[int(int(loc1)/res)]
				norm2 = kr2[int(int(loc2)/res)]
				if not np.isnan(norm1) and not np.isnan(norm2):
					out.write("\t".join((chrom1, loc1, str(int(loc1) + res), chrom2, loc2, str(int(loc2) + res), str(float(line[2])/(norm1 * norm2)))) + "\n")
		out.close()
	raw.close()

def normalize_inter(hic_id, res, chrom_a, chrom_b):
	res_string = tools.get_res_string(res)

	if get_chrom_num(chrom_a) < get_chrom_num(chrom_b):
		chrom1 = chrom_a
		chrom2 = chrom_b
	else:
		chrom1 = chrom_b
		chrom2 = chrom_a

	rawpath = "{}/{}_resolution_interchromosomal/chr{}_chr{}/MAPQGE30/chr{}_{}_{}.RAWobserved".format(hic_id, res_string, chrom1, chrom2, chrom1, chrom2, res_string)
	krpath1 = "{}/{}_resolution_interchromosomal/chr{}_chr{}/MAPQGE30/chr{}_{}.KRnorm".format(hic_id, res_string, chrom1, chrom2, chrom1, res_string)
	if not os.path.isfile(krpath1):
		krpath1 = "{}/{}_resolution_interchromosomal/chr{}_chr{}/MAPQGE30/chr{}_{}.VCnorm".format(hic_id, res_string, chrom1, chrom2, chrom1, res_string)
	krpath2 = "{}/{}_resolution_interchromosomal/chr{}_chr{}/MAPQGE30/chr{}_{}.KRnorm".format(hic_id, res_string, chrom1, chrom2, chrom2, res_string)
	if not os.path.isfile(krpath2):
		krpath2 = "{}/{}_resolution_interchromosomal/chr{}_chr{}/MAPQGE30/chr{}_{}VCnorm".format(hic_id, res_string, chrom1, chrom2, chrom2, res_string)
	outpath = "{}_{}_{}_{}.bed".format(hic_id, chrom1, chrom2, res_string)
	chromstring1 = "chr" + chrom1
	chromstring2 = "chr" + chrom2
	normalize(chromstring1, chromstring2, rawpath, krpath1, krpath2, res, outpath)

def normalize_intra(hic_id, res, chrom):
	res_string = tools.get_res_string(res)

	rawpath = "{}/{}_resolution_intrachromosomal/chr{}/MAPQGE30/chr{}_{}.RAWobserved".format(hic_id, res_string, chrom, chrom, res_string)
	krpath = "{}/{}_resolution_intrachromosomal/chr{}/MAPQGE30/chr{}_{}.KRnorm".format(hic_id, res_string, chrom, chrom, res_string)
	if not os.path.isfile(krpath):
		krpath = "{}/{}_resolution_intrachromosomal/chr{}/MAPQGE30/chr{}_{}.VCnorm".format(hic_id, res_string, chrom, chrom, res_string)
	outpath = "{}_{}_{}.bed".format(hic_id, chrom, res_string)
	chromstring = "chr" + chrom
	normalize(chromstring, chromstring, rawpath, krpath, None, res, outpath)

def main():
	parser = argparse.ArgumentParser(description="Normalize Hi-C files using Knight-Ruiz method.")
	parser.add_argument("hic_id", help="e.g. GM12878_combined")
	parser.add_argument("res", type=int, help="resolution (bp)")
	parser.add_argument("chrom1", help="first chromosome (e.g. 1)")
	parser.add_argument("--chrom2", help="second chromosome (e.g. 2)")
	args = parser.parse_args()

	if args.chrom2 is None:
		normalize_intra(args.hic_id, args.res, args.chrom1)
	else:
		normalize_inter(args.hic_id, args.res, args.chrom1, args.chrom2)

if __name__ == "__main__":
	main()
