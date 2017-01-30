set -e

RES=$1

if [ ! -e ../data ]
	then
		mkdir ../data
fi

cd ../data

if [ ! -e GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz
fi

DIR=$(($RES/1000))"kb_resolution_interchromosomal"

if [ ! -e "GM12878_combined_interchromosomal/"$DIR ]
	then
		tar xzf GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz $DIR
fi

for CHROM1 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	for CHROM2 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
	do
  		python ../figures/normalize.py GM12878_combined $RES $CHROM1 --inter -chrom2 $CHROM2
	done
done
