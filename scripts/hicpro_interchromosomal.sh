BED_PATH=$1
MAT_PATH=$2
RES_KB=$3
RES=$((RES_KB*1000))

BEDPE_PATH=all.bed

python hicpro_to_bedpe.py $BED_PATH $MAT_PATH $BEDPE_PATH

PREFIX=${MAT_PATH%.matrix}

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
NUM_CHROMS=${#CHROMS[@]}

#interchromosomal
for i in `seq $((NUM_CHROMS-1))`
do
    CHROM1=${CHROMS[$i]}

    for j in `seq 0 $((i-1))`
    do 
        CHROM2=${CHROMS[$j]}
        cat $BEDPE_PATH | awk -v chrom1=chr$CHROM1 -v chrom2=chr$CHROM2 '($1 == chrom1 && $4 == chrom2) || ($4 == chrom1 && $1 == chrom2) {print}' > ${PREFIX}_${CHROM2}_${CHROM1}_$RES_KB.bed
    done

done  

#intrachromosomal
for i in `seq 0 $((NUM_CHROMS-1))`
do
    CHROM=${CHROMS[$i]}
    cat $BEDPE_PATH | awk -v chrom=chr$CHROM '$1 == chrom && $4 == chrom {print}' > ${PREFIX}_${CHROM}_$RES_KB.bed
done 

python ../minimds_inter.py $PREFIX $RES $RES 
