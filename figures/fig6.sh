set -e

#parameters
RES=10000
RES_KB=$(($RES/1000))
DOMAIN_SIZE_PARAMETER=0.01
MIN_DOMAIN_SIZE=0.01

#get data
bash get_gm12878.sh $RES
bash get_gm12878.sh 100000

#install MOGEN
bash install_mogen.sh

#results files
MINI_OUT="minimds_"$RES_KB"kb_output.txt"
MMDS_OUT="mmds_"$RES_KB"kb_output.txt"
CMDS_OUT="cmds_"$RES_KB"kb_output.txt"
MOGEN_OUT="mogen_"$RES_KB"kb_output.txt"

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


#run MDS
for CHROM in X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
do
	BEDPATH="data/GM12878_combined_"$CHROM"_"$RES_KB"kb.bed"
	time python ../minimds.py -l "data/GM12878_combined_"$CHROM"_100kb.bed" -p $DOMAIN_SIZE_PARAMETER -m $MIN_DOMAIN_SIZE $BEDPATH >> $MINI_OUT
	time python ../minimds.py $BEDPATH >> $MMDS_OUT
	time python ../minimds.py --classical $BEDPATH >> $CMDS_OUT
	python bed_to_list.py $BEDPATH "MOGEN/examples/hiC/input/GM12878_combined_"$CHROM"_"$RES_KB"kb.tsv"
	time java -jar MOGEN/examples/hiC/3DGenerator.jar "parameters_chr"$CHROM"_"$RES_KB"kb.txt" >> $MOGEN_OUT
done

python get_chrom_sizes.py $RES_KB

#process output to get times
cat $MINI_OUT | awk '$1 == "real" {print $2}' > "minimds_"$RES_KB"kb_times.txt"
cat $MMDS_OUT | awk '$1 == "real" {print $2}' > "mmds_"$RES_KB"kb_times.txt"
cat $CMDS_OUT | awk '$1 == "real" {print $2}' > "cmds_"$RES_KB"kb_times.txt"
cat $MOGEN_OUT | awk '$1 == "real" {print $2}' > "mogen_"$RES_KB"kb_times.txt"

#plot
python fig6.py
