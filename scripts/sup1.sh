set -e

$TIME=/usr/bin/time

bash get_gm12878.sh 1000000 22
bash get_gm12878.sh 100000 22

BEDPATH=hic_data/GM12878_combined_22_100kb.bed

#Chromosome3D

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_100kb.txt

if [ ! -e $INPUT_PATH ]
	then
		python chromosome3d_input.py $BEDPATH $INPUT_PATH
fi

#run
$TIME -f %e -o chromosome3d_chr22_100kb_time.txt perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_100kb -m 1

#mMDS

$TIME -f %e -o mmds_chr22_100kb_time.txt python ../minimds.py $BEDPATH

#cMDS

$TIME -f %e -o cmds_chr22_100kb_time.txt python ../minimds.py --classical $BEDPATH

#miniMDS

$TIME -f %e -o minimds_chr22_100kb_time.txt python ../minimds.py -l hic_data/GM12878_combined_22_1mb.bed -p 0.01 -m 0.01 $BEDPATH

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
$TIME -f %e -o mogen_chr22_100kb_time.txt java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_100kb.txt

#HSA

#install
bash install_hsa.sh

#create input
INPUT_PATH=hsa/GM12878_combined_22_100kb.tsv

if [ ! -e $INPUT_PATH ]
	then
		python hsa_input.py $BEDPATH $INPUT_PATH
fi

cd hsa

#run
$TIME -f %e -o hsa_chr22_100kb_time.txt Rscript myR.R GM12878_combined_22_100kb.tsv 0 GM12878_combined_22_100kb_coords.tsv 1

cd ..

#ChromSDE

#install
#bash install_chromsde.sh

#create input
#CONTACTS_PATH=ChromSDE/chr22_100kb_contacts.dat
#IDS_PATH=ChromSDE/chr22_100kb_ids.dat

#if [ ! -e $CONTACTS_PATH ] || [ ! -e $IDS_PATH ]
#	then
#		python chromsde_input.py $BEDPATH $CONTACTS_PATH $IDS_PATH
#fi

#cd ChromSDE

#run
#$TIME -f %e -o chromsde_chr22_100kb_time.txt matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde_100kb(22)')"

#cd ..

python get_chrom_sizes.py 10

python sup1.py
