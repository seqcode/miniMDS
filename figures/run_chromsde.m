function output = run_chromsde(chrom)	
if class(chrom) == 'double'
	chrom = num2str(chrom);
end
test = sparse(importdata(strcat('chr', chrom, '_10kb_contacts.dat')));
binAnno = importdata(strcat('chr', chrom, '_10kb_ids.dat'));  
ChromSDE(binAnno, test, 1)
output = 0
exit
end

