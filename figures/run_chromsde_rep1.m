function output = run_chromsde(chrom)	
if class(chrom) == 'double'
	chrom = num2str(chrom);
end
contacts_rep1 = sparse(importdata(strcat('chr', chrom, '_10kb_contacts.dat')));
ids = importdata(strcat('chr', chrom, '_10kb_ids.dat'));  
ChromSDE(ids, contacts, 1)
output = 0
exit
end

