set -e

TIME=/usr/bin/time

BEDPATH=hic_data/GM12878_combined_22_100kb.bed

#Chromosome3D

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_100kb.txt

if [ ! -e $INPUT_PATH ]
	then
		python chromosome3d_input.py $BEDPATH $INPUT_PATH
fi

#run
$TIME -f "%M" -o chromosome3d_chr22_100kb_memory.txt perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_100kb -m 1

#mMDS

$TIME -f "%M" -o mmds_chr22_100kb_memory.txt python ../minimds.py $BEDPATH

#cMDS

$TIME -f "%M" -o cmds_chr22_100kb_memory.txt python ../minimds.py --classical $BEDPATH

#miniMDS

$TIME -f "%M" -o minimds_chr22_100kb_memory.txt python ../minimds.py -l hic_data/GM12878_combined_22_100kb.bed $BEDPATH

#MOGEN

#install
bash install_mogen.sh

#create input
INPUT_PATH=MOGEN/examples/hiC/input/GM12878_combined_22_100kb.tsv

if [ ! -e $INPUT_PATH ]
	then
		python mogen_input.py $BEDPATH $INPUT_PATH
fi

#run
$TIME -f "%M" -o mogen_chr22_100kb_memory.txt java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_100kb.txt

#ChromSDE

#install
bash install_chromsde.sh

#create input
CONTACTS_PATH=ChromSDE/chr22_100kb_contacts.dat
IDS_PATH=ChromSDE/chr22_100kb_ids.dat

if [ ! -e $CONTACTS_PATH ] || [ ! -e $IDS_PATH ]
	then
		python chromsde_input.py $BEDPATH $CONTACTS_PATH $IDS_PATH
fi

cd ChromSDE

#run
$TIME -f "%M" -o chromsde_chr22_100kb_memory.txt matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde('chr22_100kb_contacts.dat', 'chr22_100kb_ids.dat')')"

cd ..

python sup2.py
