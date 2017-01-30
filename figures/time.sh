set -e

#parameters
RES=$1
RES_KB=$(($RES/1000))
DOMAIN_SIZE_PARAMETER=0.05
MIN_DOMAIN_SIZE=0.01

#get data
bash get_gm12878.sh $RES
bash get_gm12878.sh 100000

#results files
MINI_OUT="minimds_"$RES_KB"kb_output.txt"
MMDS_OUT="mmds_"$RES_KB"kb_output.txt"
CMDS_OUT="cmds_"$RES_KB"kb_output.txt"

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


#run algorithms
for CHROM in X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
do
	INPATH="data/GM12878_combined_"$CHROM"_"$RES_KB"kb.bed"
	time python ../minimds.py -l "data/GM12878_combined_"$CHROM"_100kb.bed" -p $DOMAIN_SIZE_PARAMETER -m $MIN_DOMAIN_SIZE -o "data/GM12878_combined_"$CHROM"_"$RES_KB"kb_cluster.tsv" $INPATH >> $MINI_OUT
	time python ../minimds.py $INPATH >> $MMDS_OUT
	time python ../minimds.py --classical $INPATH >> $CMDS_OUT
done

python get_chrom_sizes.py $RES_KB

#process output to get times
cat $MINI_OUT | awk '$1 == "real" {print $2}' > "minimds_"$RES_KB"kb_times.txt"
cat $MMDS_OUT | awk '$1 == "real" {print $2}' > "mmds_"$RES_KB"kb_times.txt"
cat $CMDS_OUT | awk '$1 == "real" {print $2}' > "cmds_"$RES_KB"kb_times.txt"

#plot
python plot_time.py
