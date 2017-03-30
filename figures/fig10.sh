RES=$1

bash get_gm12878.sh $RES
bash get_gm12878.sh 100000
bash get_gm12878_inter.sh 1000000

time python ../minimds_inter.py -l 100000 -p 0.01 -m 0.01 -o data/GM12878_combined data/GM12878_combined_interchromosomal data/GM12878_combined 1000000 $RES
python plot_inter.py
