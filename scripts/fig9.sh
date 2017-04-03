set -e

bash install_mogen.sh

bash get_gm12878.sh 100000 22
bash get_gm12878.sh 10000 22

BEDPATH=hic_data/GM12878_combined_22_10kb.bed

python ../minimds.py -o hic_data/GM12878_combined_22_10kb_mmds_coords.tsv $BEDPATH
python ../minimds.py --classical -o hic_data/GM12878_combined_22_10kb_cmds_coords.tsv $BEDPATH
python ../minimds.py -l hic_data/GM12878_combined_22_100kb.bed -o hic_data/GM12878_combined_22_10kb_minimds_coords.tsv $BEDPATH

INPUT_PATH=MOGEN/examples/hiC/input/GM12878_combined_22_10kb.tsv
if [ ! -e $INPUT_PATH ]
	then
		python mogen_input.py $BEDPATH $INPUT_PATH
fi

java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_10kb.txt

#process MOGEN output
REP_NUM=1
for f in MOGEN/examples/hiC/output/GM12878_combined_22_10kb_*.pdb
do
	cat $f | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > "MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep"$REP_NUM"_coords.tsv"
	REP_NUM=$(($REP_NUM+1))
done

python fig9.py
