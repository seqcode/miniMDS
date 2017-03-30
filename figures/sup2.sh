set -e

TIME=/usr/bin/time

BEDPATH=hic_data/GM12878_combined_22_100kb.bed

#Chromosome3D

#install
bash install_chromosome3d.sh

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_100kb.txt
python chromosome3d_input.py $BEDPATH $INPUT_PATH

#run
$TIME -f "%M" -o chromosome3d_chr22_100kb_memory.txt perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output/chr22_100kb -m 1

#mMDS

$TIME -f "%M" -o mmds_chr22_100kb_memory.txt python ../minimds.py $BEDPATH

#cMDS

$TIME -f "%M" -o cmds_chr22_100kb_memory.txt python ../minimds.py --classical $BEDPATH

#miniMDS

$TIME -f "%M" -o minimds_chr22_100kb_memory.txt python ../minimds.py -l data/GM12878_combined_22_100kb.bed $BEDPATH

#MOGEN

#install
bash install_mogen.sh

#create input
python mogen_input.py $BEDPATH MOGEN/examples/hiC/input/GM12878_combined_22_100kb.tsv

#run
$TIME -f "%M" -o mogen_chr22_100kb_memory.txt java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_100kb.txt

#ChromSDE

#install
bash install_chromsde.sh

#create input
python chromsde_input.py $BEDPATH ChromSDE/chr22_100kb_contacts.dat ChromSDE/chr22_100kb_ids.dat

cd ChromSDE

#run
$TIME -f "%M" -o chromsde_chr22_100kb_memory.txt matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde('chr22_100kb_contacts.dat', 'chr22_100kb_ids.dat')')"

cd ..

python sup2.py
