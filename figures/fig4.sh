#results files
CHROMOSOME3D_OUT="chromosome3d_chr22_10kb_output.txt"
MOGEN_OUT="mogen_chr22_10kb_output.txt"
MINI_OUT="minimds_chr22_10kb_output.txt"
MMDS_OUT="mmds_chr22_10kb_output.txt"
CMDS_OUT="cmds_chr22_10kb_output.txt"
CHROMSDE_OUT="chromsde_chr22_10kb_output.txt"

bash get_gm12878.sh 100000
bash get_gm12878.sh 10000

BEDPATH=hic_data/GM12878_combined_22_10kb.bed

#Chromosome3D

#install
bash install_chromosome3d.sh

#create input
INPUT_PATH=Chromosome3D/input/GM12878_combined_22_10kb.txt
python chromosome3d_input.py $BEDPATH $INPUT_PATH

#run
time perl Chromosome3D/chromosome3D.pl -i $INPUT_PATH -o Chromosome3D/output/chr22_10kb -m 1 > $CHROMOSOME3D_OUT

#mMDS

time python ../minimds.py $BEDPATH > $MMDS_OUT

#cMDS

time python ../minimds.py --classical $BEDPATH > $CMDS_OUT

#miniMDS

time python ../minimds.py -l data/GM12878_combined_22_100kb.bed -p 0.01 -m 0.01 $BEDPATH > $MINI_OUT

#MOGEN

#install
bash install_mogen.sh

#create input
python mogen_input.py $BEDPATH MOGEN/examples/hiC/input/GM12878_combined_22_10kb.tsv

#run
time java -jar MOGEN/examples/hiC/3DGenerator.jar parameters_chr22_10kb.txt > $MOGEN_OUT

#ChromSDE

#install
bash install_chromsde.sh

cd ChromSDE

#create input
python chromsde_input.py $BEDPATH chr22_10kb_contacts.dat chr22_10kb_ids.dat

#run
time matlab -nodisplay -nosplash -nodesktop -r "run('run_chromsde(22)')" > "../"$CHROMSDE_OUT

cd ..

python get_chrom_sizes.py 10

#process output to get times
TIMES=chr22_10kb_times.txt
cat $CHROMOSOME3D_OUT | awk '$1 == "real" {print $2}' > $TIMES
cat $MINI_OUT | awk '$1 == "real" {print $2}' >> $TIMES
cat $MMDS_OUT | awk '$1 == "real" {print $2}' >> $TIMES
cat $CMDS_OUT | awk '$1 == "real" {print $2}' >> $TIMES
cat $MOGEN_OUT | awk '$1 == "real" {print $2}' >> $TIMES
cat $CHROMSDE_OUT | awk '$1 == "real" {print $2}' >> $TIMES

python fig4.py
