set -e

#parameters
RES=10000
RES_KB=$(($RES/1000))
TIME=/usr/bin/time
DOMAIN_SIZE_PARAMETER=0.01
MIN_DOMAIN_SIZE=0.01

#install MOGEN
bash install_mogen.sh

#get data
bash get_gm12878.sh $RES
bash get_gm12878.sh 100000

#results files
MINI_OUT="minimds_"$RES_KB"kb_memory.txt"
MMDS_OUT="mmds_"$RES_KB"kb_memory.txt"
CMDS_OUT="cmds_"$RES_KB"kb_memory.txt"
MOGEN_OUT="mogen_"$RES_KB"kb_memory.txt"

#reset
if [ -e $MINI_OUT ]
	then
		rm $MINI_OUT
fi

if [ -e $MMDS_OUT ]
	then
		rm $MMDS_OUT
fi

if [ -e $CMDS_OUT ]
	then
		rm $CMDS_OUT
fi

if [ -e $MOGEN_OUT ]
	then
		rm $MOGEN_OUT
fi

#run algorithms
for CHROM in X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
do
	INPATH="data/GM12878_combined_"$CHROM"_"$RES_KB"kb.bed"
       	$TIME -f "%M" -o $CMDS_OUT -a python ../minimds.py --classical $INPATH
	$TIME -f "%M" -o $MMDS_OUT -a python ../minimds.py $INPATH
	$TIME -f "%M" -o $MINI_OUT -a python ../minimds.py -l "data/GM12878_combined_"$CHROM"_100kb.bed" -p $DOMAIN_SIZE_PARAMETER -m $MIN_DOMAIN_SIZE $INPATH 
	$TIME -f "%M" -o $MINI_OUT -a java -jar MOGEN/examples/hiC/3DGenerator.jar "parameters_chr"$CHROM"_"$RES_KB"kb.txt"
done

python get_chrom_sizes.py $RES_KB

#plot
python plot_memory.py
