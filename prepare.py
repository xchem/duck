import sys
import openmoltools as moltools
import openbabel
from rdkit.Chem import AllChem
from rdkit import Chem
# Set the names
prefix = sys.argv[1]
mol_file = prefix + ".mol"
# Get the charge
rd_mol = Chem.MolFromMolFile(mol_file)
net_charge = AllChem.GetFormalCharge(rd_mol)
mol2_file = prefix + ".mol2"
out_file = prefix + "_params.mol2"
frcmod_file = prefix + ".frcmod"
# Conver the molecule
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("mdl","mol2")
mol = openbabel.OBMol()
obConversion.ReadFile(mol, mol_file)
mol.AddHydrogens()
obConversion.WriteFile(mol, mol2_file)
moltools.amber.run_antechamber(prefix, mol2_file,net_charge=net_charge,gaff_mol2_filename=out_file,frcmod_filename=frcmod_file)
