import sys,os
import openbabel
from rdkit.Chem import AllChem
from rdkit import Chem

def prep_lig(mol_file,prefix):
    # Get the charge
    rd_mol = Chem.MolFromMolFile(mol_file)
    net_charge = AllChem.GetFormalCharge(rd_mol)
    mol2_file = prefix + ".mol2"
    out_file = prefix + "_params.mol2"
    frcmod_file = prefix + ".frcmod"
    # Conver the molecule
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mdl", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, mol_file)
    mol.AddHydrogens()
    obConversion.WriteFile(mol, mol2_file)
    os.system("antechamber -i " + mol2_file + " -fi mol2 -o " + out_file + " -fo mol2 -at sybyl  -c bcc")
    return [mol2_file]

if __name__ == "__main__":
    # Set the names
    prep_lig(sys.argv[2],sys.argv[1])