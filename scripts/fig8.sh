set -e

bash install_mogen.sh

bash get_gm12878.sh 100000 0
bash get_gm12878.sh 10000 0

for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	BEDPATH="hic_data/GM12878_combined_"$CHROM"_10kb.bed"
	python ../minimds.py -o "hic_data/GM12878_combined_"$CHROM"_10kb_mmds_coords.tsv" $BEDPATH
	python ../minimds.py --classical -o "hic_data/GM12878_combined_"$CHROM"_10kb_cmds_coords.tsv" $BEDPATH
	python ../minimds.py -l "hic_data/GM12878_combined_"$CHROM"_100kb.bed" -o "hic_data/GM12878_combined_"$CHROM"_10kb_minimds_coords.tsv" -p 0.001 $BEDPATH

	INPUT_PATH="MOGEN/examples/hiC/input/GM12878_combined_"$CHROM"_"$RES_KB"kb.tsv" 
	if [ ! -e $INPUT_PATH ]
		then
			python mogen_input.py $BEDPATH $INPUT_PATH
	fi
	java -jar MOGEN/examples/hiC/3DGenerator.jar "parameters_chr"$CHROM"_10kb.txt"

	#process MOGEN output
	REP_NUM=1
	for f in "MOGEN/examples/hiC/output/GM12878_combined_"$CHROM"_10kb_"*".pdb"
	do
		cat $f | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > "MOGEN/examples/hiC/output/GM12878_combined_"$CHROM"_10kb_rep"$REP_NUM"_coords.tsv"
		REP_NUM=$(($REP_NUM+1))
	done
done

python fig8.py
