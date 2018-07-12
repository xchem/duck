from simtk.openmm import app
from rdkit import Chem
import parmed
import sys
import pickle
from duck.utils.gen_system import generateSMIRNOFFStructureRDK


def write_pickle(out_file,out_pmd):
    pickle_out = open(out_file, "wb")
    pickle.dump([out_pmd], pickle_out)
    pickle_out.close()

def prepare_system(ligand_file, protein_file, forcefield_str='amber99sb.xml', water_ff='amber99_obc.xml'):
    print("Preparing ligand")
    lig_rdk = Chem.MolFromMol2File(ligand_file, sanitize=False)
    lig_rdk.SetProp('_Name', 'LIG')
    ligand_pmd = generateSMIRNOFFStructureRDK(lig_rdk)
    print("Fixing protein")
    protein = parmed.load_file(protein_file)
    protein.write_pdb("fixed.pdb")
    print("loading system")
    protein = parmed.load_file("fixed.pdb")["!(:HOH,NA,CL)"]  # remove ions and water
    forcefield = app.ForceField(forcefield_str,water_ff)
    protein_system = forcefield.createSystem(protein.topology)
    protein_pmd = parmed.openmm.load_topology(protein.topology, protein_system, protein.positions)
    protein_pmd.save("protein_prepared.pdb", overwrite=True)
    prot_lig_pmd = protein_pmd + ligand_pmd
    prot_lig_pmd.save("complex.pdb", overwrite=True)
    print("Parametrizing solvent done")
    print("merge structures")
    combined_pmd = protein_pmd + ligand_pmd
    combined_pmd.save("system_complex.prmtop", overwrite=True)
    print("merge done")
    combined_pmd.box_vectors = complex.box_vectors
    print("writing pickles")
    complex_out = "complex_system.pickle"
    protein_out = "protein_system.pickle"
    ligand_out = "ligand_system.pickle"
    write_pickle(complex_out,complex)
    write_pickle(ligand_out,ligand_pmd)
    write_pickle(protein_out,protein_system)
    return [complex_out, protein_out, ligand_out]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE : python prepare_mmgbsa.py protein.pdb ligand.mol2\n")
    ligand_file = sys.argv[2]
    protein_file = sys.argv[1]
    prepare_system(ligand_file, protein_file)
