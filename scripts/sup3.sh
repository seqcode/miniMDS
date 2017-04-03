set -e

bash get_gm12878.sh 1000000 22
bash get_gm12878.sh 100000 22

BEDPATH=hic_data/GM12878_combined_22_100kb.bed

#Chromosome3D

#install

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_100kb.txt

if [ ! -e $INPUT_PATH ]
	then
		python chromosome3d_input.py $BEDPATH $INPUT_PATH
fi

#run
perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_100kb -m 1

#process output
cat Chromosome3D/output_models/chr22_100kb/GM12878_combined_22_100kb_1.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > Chromosome3D/output_models/chr22_100kb/chr22_100kb_coords.tsv

#mMDS

python ../minimds.py $BEDPATH

#cMDS

python ../minimds.py --classical $BEDPATH

#miniMDS

python ../minimds.py -l hic_data/GM12878_combined_22_1mb.bed -p 0.01 -m 0.01 $BEDPATH

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
java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_100kb.txt

#process output
REP_NUM=1
for f in MOGEN/examples/hiC/output/GM12878_combined_22_100kb_*.pdb
do
	cat $f | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > "MOGEN/examples/hiC/output/GM12878_combined_22_100kb_rep"$REP_NUM"_coords.tsv"
	REP_NUM=$(($REP_NUM+1))
done

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
Rscript myR.R GM12878_combined_22_100kb.tsv 0 GM12878_combined_22_100kb_coords.tsv 1

cd ..

#ChromSDE

#install
bash install_chromsde.sh

#create input
CONTACTS_PATH=ChromSDE/chr22_100kb_contacts.dat
IDS_PATH=ChromSDE/chr22_100kb_ids.dat

python chromsde_input.py $BEDPATH $CONTACTS_PATH $IDS_PATH

cd ChromSDE

#run
matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde_100kb(22)')"

#process output
cat contacts_100kb.pos.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > GM12878_combined_22_100kb_coords.tsv

cd ..

python sup3.py
