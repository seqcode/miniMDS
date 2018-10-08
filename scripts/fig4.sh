set -e

TIME=/usr/bin/time

bash get_gm12878.sh 100000 22
bash get_gm12878.sh 10000 22

BEDPATH=hic_data/GM12878_combined_22_10kb.bed

#Chromosome3D

#create input
#INPUT_PATH=Chromosome3D/input/GM12878_combined_22_10kb.txt

#if [ ! -e $INPUT_PATH]
#	then
#		python chromosome3d_input.py $BEDPATH $INPUT_PATH
#fi

#run
#$TIME -o chromosome3d_chr22_10kb_time.txt -f %e perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_10kb -m 1

#mMDS

$TIME -o mmds_chr22_10kb_time.txt -f %e python ../minimds.py $BEDPATH

#cMDS

$TIME -o cmds_chr22_10kb_time.txt -f %e python ../minimds.py --classical $BEDPATH

#miniMDS

$TIME -o minimds_chr22_10kb_time.txt -f %e python ../minimds.py --partitioned $BEDPATH

#MOGEN

#install
bash install_mogen.sh

#create input
INPUT_PATH=MOGEN/examples/hiC/input/GM12878_combined_22_10kb.tsv
if [ ! -e $INPUT_PATH ]
	then
		python mogen_input.py $BEDPATH $INPUT_PATH
fi

#run
$TIME -o mogen_chr22_10kb_time.txt -f %e java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_10kb.txt

#ChromSDE

#install
#bash install_chromsde.sh

#create input
#CONTACTS_PATH=ChromSDE/chr22_10kb_contacts.dat
#IDS_PATH=ChromSDE/chr22_10kb_ids.dat

#if [ ! -e $CONTACTS_PATH ] || [ ! -e $IDS_PATH ]
#	then
#		python chromsde_input.py $BEDPATH $CONTACTS_PATH $IDS_PATH
#fi

#cd ChromSDE

#run
#$TIME -o chromsde_chr22_10kb_time.txt -f %e matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde('chr22_10kb_contacts.dat', 'chr22_10kb_ids.dat')')" > "../"$CHROMSDE_OUT

#cd ..

python fig4.py
