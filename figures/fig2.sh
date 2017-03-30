if [ ! -e mESC_chr6.tsv ]
	then
		wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35156/suppl/GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt.gz
		gunzip GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt.gz
		cat GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt | awk ' $2 == "chr6" && $5 == "chr6" {print $3"\t"$6} ' > mESC_chr6.tsv
		rm GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt
fi
python fig2.py
