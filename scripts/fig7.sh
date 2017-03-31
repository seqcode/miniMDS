set -e

#parameters
RES=10000
RES_KB=$(($RES/1000))
TIME=/usr/bin/time
DOMAIN_SIZE_PARAMETER=0.01
MIN_DOMAIN_SIZE=0.01

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

#install MOGEN
bash install_mogen.sh

#get data
bash get_gm12878.sh $RES 0
bash get_gm12878.sh 100000 0

#run algorithms
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	BEDPATH="data/GM12878_combined_"$CHROM"_"$RES_KB"kb.bed"
       	$TIME -f "%M" -o $CMDS_OUT -a python ../minimds.py --classical $BEDPATH
	$TIME -f "%M" -o $MMDS_OUT -a python ../minimds.py $BEDPATH
	$TIME -f "%M" -o $MINI_OUT -a python ../minimds.py -l "data/GM12878_combined_"$CHROM"_100kb.bed" -p $DOMAIN_SIZE_PARAMETER -m $MIN_DOMAIN_SIZE $BEDPATH

	INPUT_PATH="MOGEN/examples/hiC/input/GM12878_combined_"$CHROM"_"$RES_KB"kb.tsv" 
	if [ ! -e $INPUT_PATH ]
		then
			python mogen_input.py $BEDPATH $INPUT_PATH
	fi
	$TIME -f "%M" -o $MINI_OUT -a java -jar MOGEN/examples/hiC/3DGenerator.jar "parameters_chr"$CHROM"_"$RES_KB"kb.txt"
done

if [ ! -e "chrom_sizes_"$RES_KB"kb.txt" ]
	then
		python get_chrom_sizes.py $RES_KB
fi

#plot
python fig7.py
