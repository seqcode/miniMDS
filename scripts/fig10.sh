set -e

#./get_gm12878.sh 10000 0 
#./get_gm12878_inter.sh 1000000

python ../minimds_inter.py --partitioned -l 10 -p 0.01 -m 0.01 -o hic_data/GM12878_combined hic_data/GM12878_combined 1000000 10000
python fig10.py
