function output = run_chromsde(chrom)	
if class(chrom) == 'double'
	chrom = num2str(chrom);
end
contacts_100kb = sparse(importdata(strcat('chr', chrom, '_100kb_contacts.dat')));
ids = importdata(strcat('chr', chrom, '_100kb_ids.dat'));  
ChromSDE(ids, contacts, 1)
output = 0
exit
end

