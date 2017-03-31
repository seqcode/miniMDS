if [ ! -e ChromSDE ]
	then
		wget http://biogpu.ddns.comp.nus.edu.sg/~chipseq/ChromSDE/ChromSDE_program2.2.zip
		unzip ChromSDE_program2.2.zip
		rm ChromSDE_program2.2.zip
		mv program ChromSDE
		mv *.m ChromSDE
fi
