from simtk.openmm import app
from rdkit import Chem
from simtk import unit
import parmed, pkg_resources
from simtk import openmm
from pdbfixer import PDBFixer  # for solvating
import sys
import pickle
from duck.utils.gen_system import generateSMIRNOFFStructureRDK


def find_box_size(input_file="complex.pdb",add_factor=20):
    complex = parmed.load_file(input_file)
    x_coords = [x[0] for x in complex.positions]
    y_coords = [x[1] for x in complex.positions]
    z_coords = [x[2] for x in complex.positions]
    x_size = abs(max(x_coords)-min(x_coords))
    y_size = abs(max(y_coords)-min(y_coords))
    z_size = abs(max(z_coords)-min(z_coords))
    val_in_ang = max(x_size,y_size,z_size) + add_factor*unit.angstrom
    return int(val_in_ang.value_in_unit(unit.angstrom))+1


def prepare_system(ligand_file, protein_file, forcefield_str='amber99sb.xml'):
    print("Preparing ligand")
    lig_rdk = Chem.MolFromMol2File(ligand_file, sanitize=False)
    lig_rdk.SetProp('_Name', 'LIG')
    ligand_pmd = generateSMIRNOFFStructureRDK(lig_rdk)
    print("Fixing protein")
    protein = parmed.load_file(protein_file)
    protein.write_pdb("fixed.pdb")
    print("loading system")
    protein = parmed.load_file("fixed.pdb")["!(:HOH,NA,CL)"]  # remove ions and water
    forcefield = app.ForceField(forcefield_str)
    protein_system = forcefield.createSystem(protein.topology)
    protein_pmd = parmed.openmm.load_topology(protein.topology, protein_system, protein.positions)
    protein_pmd.save("protein_prepared.pdb", overwrite=True)
    prot_lig_pmd = protein_pmd + ligand_pmd
    prot_lig_pmd.save("complex.pdb", overwrite=True)
    print("Solvation")
    fixer = PDBFixer("complex.pdb")
    # 0.1 in Vec3 because box_size is in angstrom and fixer uses nanometer
    # scaling factor to somehow ensure no interaction with periodic image
    scaling_factor = 1.0
    box_size = find_box_size("complex.pdb")
    fixer.addSolvent(scaling_factor * box_size * openmm.Vec3(0.1, 0.1, 0.1), positiveIon='Na+', negativeIon='Cl-',
                     ionicStrength=0.1 * unit.molar)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open('complex_solvated.pdb', 'w'))
    print("Solvation done")
    print("Parametrizing ions")
    complex = parmed.load_file('./complex_solvated.pdb')
    ions = complex["(:NA,CL)"]
    forcefield = app.ForceField(forcefield_str)
    ions_system = forcefield.createSystem(ions.topology)
    ions_pmd = parmed.openmm.load_topology(ions.topology, ions_system, ions.positions)
    print("Parametrizing ions done")
    prot_lig_ion_pmd = protein_pmd + ligand_pmd + ions_pmd
    print("Parametrizing solvent")
    solvent = complex["(:HOH)"]
    num_solvent = len(solvent.residues)
    prm_top_water_path = pkg_resources.resource_filename('duck', "parameters/waters/water.prmtop")
    solvent_pmd = parmed.load_file(prm_top_water_path)
    solvent_pmd *= num_solvent
    solvent_pmd.positions = solvent.positions
    print("Parametrizing solvent done")
    print("merge structures")
    combined_pmd = protein_pmd + ligand_pmd + ions_pmd + solvent_pmd
    combined_pmd.save("system_complex.prmtop", overwrite=True)
    print("merge done")
    combined_pmd.box_vectors = complex.box_vectors
    print("writing pickle")
    complex = "./complex_system.pickle"
    pickle_out = open(complex, "wb")
    pickle.dump([combined_pmd], pickle_out)
    pickle_out.close()
    return [complex]


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("USAGE : python parametrize.py protein.pdb ligand.mol2\n")
    ligand_file = sys.argv[2]
    protein_file = sys.argv[1]
    prepare_system(ligand_file, protein_file)
