BED_PATH=$1
MAT_PATH=$2
RES=$3

BEDPE_PATH=all.bed

python hicpro_to_bedpe.py $BED_PATH $MAT_PATH $BEDPE_PATH

PREFIX=${MAT_PATH%.matrix}