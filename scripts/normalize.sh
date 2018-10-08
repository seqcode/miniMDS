set -e

for CHROM in 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20
do
	echo $CHROM
	python normalize.py /data/drive1/hic_data/K562 25000 $CHROM
done
