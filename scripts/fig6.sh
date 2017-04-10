set -e

TIME=/usr/bin/time

#parameters
DOMAIN_SIZE_PARAMETER=0.01
MIN_DOMAIN_SIZE=0.01

#results files
MINI_OUT=minimds_10kb_times.txt
MMDS_OUT=mmds_10kb_times.txt
CMDS_OUT=cmds_10kb_times.txt
MOGEN_OUT=mogen_10kb_times.txt

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

#get data
bash get_gm12878.sh 10000 0
bash get_gm12878.sh 100000 0

#install MOGEN
bash install_mogen.sh

#run MDS
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	BEDPATH="hic_data/GM12878_combined_"$CHROM"_10kb.bed"
	$TIME -o $MINI_OUT -a -f %e python ../minimds.py -l "hic_data/GM12878_combined_"$CHROM"_100kb.bed" -p $DOMAIN_SIZE_PARAMETER -m $MIN_DOMAIN_SIZE $BEDPATH
	$TIME -o $MMDS_OUT -a -f %e python ../minimds.py $BEDPATH
	$TIME -o $CMDS_OUT -a -f %e python ../minimds.py --classical $BEDPATH


	INPUT_PATH="MOGEN/examples/hiC/input/GM12878_combined_"$CHROM"_10kb.tsv"
	if [ ! -e $INPUT_PATH ]
		then
			python mogen_input.py $BEDPATH $INPUT_PATH
	fi
	$TIME -o $MOGEN_OUT -a f %e java -jar MOGEN/examples/hiC/3DGenerator.jar "parameters_chr"$CHROM"_10kb.txt"
done

if [ ! -e chrom_sizes_10kb.txt ] || [ ! -s chrom_size_10kb.txt ]
	then
		python get_chrom_sizes.py 10
fi

#plot
python fig6.py
