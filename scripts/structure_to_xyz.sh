f=$1
PREFIX=${f%.tsv}

cat $f | awk 'NF == 4 && $2 != "nan" {print "C\t"$2"\t"$3"\t"$4}' > $PREFIX.xyz