from Bio import PDB
import sys
sys.path.append("..")
from data_tools import structure_from_file

pdb_path = sys.argv[1]
scale = float(sys.argv[2])
struct_paths = sys.argv[3:]

# Create a new PDB structure
structure = PDB.Structure.Structure(0)

# Add a model to the structure
model = PDB.Model.Model(0)
structure.add(model)

atom_num = 0
all_atom_nums = []
for i, path in enumerate(struct_paths):
    # Add a chain to the model
    chain = PDB.Chain.Chain(str(i))
    model.add(chain)

    # Add a residue to the chain
    residue = PDB.Residue.Residue((0,0,0), 0, 0)
    chain.add(residue)
    
    atom_nums = []
    for x,y,z in structure_from_file(path).getCoords():
        # Add an atom to the residue
        residue.add(PDB.Atom.Atom(atom_num, (x*scale, y*scale, z*scale), 0, 0, 0, "CA", 0, "C"))
        atom_num += 1
        atom_nums.append(atom_num)
    all_atom_nums.append(atom_nums)

# Write the PDB file
io = PDB.PDBIO()
io.set_structure(structure)
io.save(pdb_path)

#add connections
with open(pdb_path, "a") as out_file:
    for atom_nums in all_atom_nums:
        for i in range(1,len(atom_nums)):
            atom_num1 = atom_nums[i-1]
            atom_num2 = atom_nums[i]
            out_file.write(f"CONECT{atom_num1:>5}{atom_num2:>5}\n")