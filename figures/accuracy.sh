set -e

RES=$1

bash get_gm12878.sh $RES
bash get_gm12878.sh 100000
python mds_accuracy.py $RES
