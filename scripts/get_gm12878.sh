set -e

RES=$1
CHROM=$2
RES_KB=$(($RES/1000))

if [ ! -e hic_data ]
	then
		mkdir hic_data
fi

cd hic_data

if [ ! -e GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Fintrachromosomal%5Fcontact%5Fmatrices%2Etar%2Egz
fi

if [ $CHROM -eq 0 ]
	then
		DIR=$RES_KB"kb_resolution_intrachromosomal"
	else			
		DIR=$RES_KB"kb_resolution_intrachromosomal/chr"$CHROM
fi

if [ ! -e "GM12878_combined/"$DIR ]
	then
		tar xzf GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz "GM12878_combined/"$DIR
fi

if [ $CHROM -eq 0 ]
	then
		for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
		do
			if [ ! -e "GM12878_combined_"$CHROM"_"$RES_KB"kb.bed" ]
				then
	  				python ../normalize.py GM12878_combined $RES $CHROM
			fi
		done

else
	if [ ! -e "GM12878_combined_"$CHROM"_"$RES_KB"kb.bed" ]
		then
	  		python ../normalize.py GM12878_combined $RES $CHROM
	fi

fi

cd ..
