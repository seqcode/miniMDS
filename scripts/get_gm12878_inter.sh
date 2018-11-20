set -e

RES=$1

mkdir -p hic_data

cd hic_data

if [ ! -e GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz ]
	then
		curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Finterchromosomal%5Fcontact%5Fmatrices%2Etar%2Egz -o GSE65325_GM12878_combined_interchromosomal_contact_matrices.tar.gz
fi

RES_KB=$(($RES/1000))

if [ $RES_KB -lt 1000 ]
	then
		RES_STRING=$RES_KB"kb"
else
	RES_STRING=$(($RES_KB/1000))"mb"
fi

DIR=$RES_STRING"_resolution_interchromosomal"

if [ ! -e GM12878_combined_interchromosomal/$DIR ]
	then
		tar xzf GSE63525_GM12878_combined_interchromosomal_contact_matrices.tar.gz GM12878_combined_interchromosomal/$DIR
fi

cd ..

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

for i in `seq 0 $((${#CHROMS[@]}-1))`
do
	CHROM1=${CHROMS[$i]}
	for j in `seq 0 $(($i-1))`
	do
		CHROM2=${CHROMS[$j]}
		python normalize.py hic_data/GM12878_combined_interchromosomal $RES $CHROM1 --chrom2 $CHROM2
		mv hic_data/GM12878_combined_interchromosomal_${CHROM2}_${CHROM1}_${RES_STRING}.bed hic_data/GM12878_combined_${CHROM2}_${CHROM1}_${RES_STRING}.bed
	done
done
