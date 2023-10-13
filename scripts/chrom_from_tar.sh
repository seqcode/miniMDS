RES=$1
CHROM=$2
TAR=$3
DATA_DIR=$4
PREFIX=$5

RES_KB=$(($RES/1000))

if [ $RES_KB -lt 1000 ]
	then
		RES_STRING=$RES_KB"kb"
else
	RES_STRING=$(($RES_KB/1000))"mb"
fi

OUT_DIR=$PREFIX/$RES_STRING"_resolution_intrachromosomal"/chr$CHROM
test ! -d $DATA_DIR/$OUT_DIR && (tar -C $DATA_DIR -xzf $TAR $OUT_DIR)

test ! -s $DATA_DIR/$PREFIX"_"$CHROM"_"$RES_KB$RES_STRING.bed && (python3 normalize.py $DATA_DIR/$PREFIX $RES $CHROM)