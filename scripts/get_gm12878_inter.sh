set -e

RES=$1

if [ ! -e hic_data ]
	then
		mkdir hic_data
fi

cd hic_data

if [ ! -e GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz
fi

RES_KB=$(($RES/1000))

if [ $RES_KB -lt 1000 ]
	then
		RES_STRING=$RES_KB"kb"
else
	RES_STRING=$(($RES_KB/1000))"mb"
fi

DIR=$RES_STRING"_resolution_interchromosomal"

if [ ! -e "GM12878_combined_interchromosomal/"$DIR ]
	then
		tar xzf GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz "GM12878_combined_interchromosomal/"$DIR
fi

for CHROM1 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	for CHROM2 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
	do
		if [ ! $CHROM1 == $CHROM2 ]
			then
  				python ../normalize.py GM12878_combined_interchromosomal $RES $CHROM1 -chrom2 $CHROM2
		fi
	done
done
