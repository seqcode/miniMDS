bash get_gm12878.sh 5000 22 
python ../minimds.py hic_data/GM12878_combined_22_5kb.bed
python ../minimds.py --partitioned hic_data/GM12878_combined_22_5kb.bed
