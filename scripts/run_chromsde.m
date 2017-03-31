function output = run_chromsde(contacts_path, ids_path)	
test = sparse(importdata(contacts_path));
binAnno = importdata(ids_path);  
ChromSDE(binAnno, test, 1)
output = 0
exit
end

