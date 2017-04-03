set -e

bash get_gm12878.sh 100000 22
bash get_gm12878.sh 10000 22

BEDPATH=hic_data/GM12878_combined_22_10kb.bed

#Chromosome3D

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_10kb.txt

if [ ! -e $INPUT_PATH ]
	then
		python chromosome3d_input.py $BEDPATH $INPUT_PATH
fi

#rep 1
perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_10kb_rep1 -m 1

#rep 2
perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output_models/chr22_10kb_rep2 -m 1

#process output
cat Chromosome3D/output_models/chr22_10kb_rep1/GM12878_combined_22_10kb_1.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > Chromosome3D/output_models/chr22_10kb_rep1/rep1_coords.tsv
cat Chromosome3D/output_models/chr22_10kb_rep2/GM12878_combined_22_10kb_1.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > Chromosome3D/output_models/chr22_10kb_rep2/rep2_coords.tsv

#mMDS

#rep 1
python ../minimds.py -o hic_data/GM12878_combined_22_10kb_mmds_rep1.tsv $BEDPATH

#rep 2
python ../minimds.py -o hic_data/GM12878_combined_22_10kb_mmds_rep2.tsv $BEDPATH

#miniMDS

#rep 1
python ../minimds.py -l hic_data/GM12878_combined_22_100kb.bed -p 0.01 -m 0.01 -o hic_data/GM12878_combined_22_10kb_minimds_rep1.tsv $BEDPATH

#rep 2
python ../minimds.py -l hic_data/GM12878_combined_22_100kb.bed -p 0.01 -m 0.01 -o hic_data/GM12878_combined_22_10kb_minimds_rep2.tsv $BEDPATH

#MOGEN

#install
bash install_mogen.sh

#create input
INPUT_PATH=MOGEN/examples/hiC/input/GM12878_combined_22_10kb.tsv
if [ ! -e $INPUT_PATH]
	then
		python mogen_input.py $BEDPATH $INPUT_PATH
fi

#rep 1
java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_10kb.txt

#rep 2
java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_10kb.txt

#process output
REP_NUM=1
for f in MOGEN/examples/hiC/output/GM12878_combined_22_10kb_*.pdb
do
	cat $f | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > "MOGEN/examples/hiC/output/GM12878_combined_22_10kb_rep"$REP_NUM"_coords.tsv"
	REP_NUM=$(($REP_NUM+1))
done

#HSA

#install
bash install_hsa.sh

#create input
INPUT_PATH=hsa/GM12878_combined_22_10kb.tsv

if [ ! -e $INPUT_PATH ]
	then
		python hsa_input.py $BEDPATH $INPUT_PATH
fi

cd hsa

#rep 1
Rscript myR.R GM12878_combined_22_10kb.tsv 0 GM12878_combined_22_10kb_rep1_coords.tsv 1

#rep 2
Rscript myR.R GM12878_combined_22_10kb.tsv 0 GM12878_combined_22_10kb_rep2_coords.tsv 1

cd ..

#ChromSDE

#install
bash install_chromsde.sh

#create input
CONTACTS_PATH=ChromSDE/chr22_10kb_contacts.dat
IDS_PATH=ChromSDE/chr22_10kb_ids.dat

if [ ! -e $CONTACTS_PATH ] || [ ! -e $IDS_PATH ]
	then
		python chromsde_input.py $BEDPATH $CONTACTS_PATH $IDS_PATH
fi

cd ChromSDE

#rep 1
matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde_rep1')"

#rep 2
matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde_rep2')"

#process output
cat contacts_rep1.pos.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > GM12878_combined_22_10kb_rep1_coords.tsv
cat contacts_rep2.pos.pdb | awk '$1 == "ATOM" {print $6"\t"$7"\t"$8}' > GM12878_combined_22_10kb_rep2_coords.tsv

cd ..

python fig5.py
