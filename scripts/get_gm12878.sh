set -e

RES=$1
CHROM=$2

DATA_DIR=hic_data
mkdir -p $DATA_DIR

PREFIX=GM12878_combined
TAR=$DATA_DIR/GSE63525_${PREFIX}_intrachromosomal_contact_matrices.tar.gz

test ! -s $TAR && (curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Fintrachromosomal%5Fcontact%5Fmatrices%2Etar%2Egz -o $TAR)

#all human chromosomes
if [ $CHROM -eq 0 ]
	then
		for CHROM in `seq 23`
		do
			bash chrom_from_tar.sh $RES $CHROM $TAR $DATA_DIR $PREFIX
		done

		bash chrom_from_tar.sh $RES X $TAR $DATA_DIR $PREFIX

#selected chromosome
else
	bash chrom_from_tar.sh $RES $CHROM $TAR $DATA_DIR $PREFIX
fi