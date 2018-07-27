import sys,os
from rdkit.Chem import AllChem
from rdkit import Chem

def prep_lig(mol_file,prefix,addHs=True):
    # Accept readily protonated molecule
    rd_mol = Chem.MolFromMolFile(mol_file,removeHs=False)
    # Add extra hydrogens
    rd_mol = AllChem.AddHs(rd_mol, addCoords=True)
    net_charge = AllChem.GetFormalCharge(rd_mol)
    # Defined the final out file
    out_file = prefix + "_params.mol2"
    frcmod_file = prefix + ".frcmod"
    os.system("antechamber -i " + mol_file + " -fi mdl -o " + out_file + " -fo mol2 -at sybyl  -c bcc -nc "+str(net_charge))
    return [out_file]

if __name__ == "__main__":
    # Set the names
    prep_lig(sys.argv[2],sys.argv[1])