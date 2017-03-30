set -e

bash get_gm12878.sh 10000 0
bash get_gm12878.sh 100000 0
bash get_gm12878_inter.sh 1000000

python ../minimds_inter.py -l 100000 -p 0.01 -m 0.01 -o hic_data/GM12878_combined data/GM12878_combined_interchromosomal hic_data/GM12878_combined 1000000 10000
python fig10.py
