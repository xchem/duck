import numpy as np
import openforcefield.utils as utils
from openforcefield.typing.engines.smirnoff import forcefield_rdk
from simtk import unit
import parmed, pkg_resources

def create_system_from_molecule_rdk(forcefield, mol, verbose=False):
	"""
	Generate a System from the given OEMol and SMIRNOFF forcefield, return the resulting System.
	Parameters
	----------
	forcefield : ForceField
	    SMIRNOFF forcefield
	mol : RDKit molecule
	    Molecule to test (must have coordinates)
	Returns
	----------
	topology : OpenMM Topology
	system : OpenMM System
	positions : initial atomic positions (OpenMM)
	"""
	# Create system
	topology = utils.generateTopologyFromRDKMol(mol)
	system = forcefield.createSystem(topology, [mol], verbose=verbose)
	# Get positions
	coordinates = mol.GetConformer().GetPositions()
	natoms = len(coordinates)
	positions = np.zeros([natoms,3], np.float32)
	for index in range(natoms):
	    (x,y,z) = coordinates[index]
	    positions[index,0] = x
	    positions[index,1] = y
	    positions[index,2] = z
	positions = unit.Quantity(positions, unit.angstroms)
	return topology, system, positions

def generateSMIRNOFFStructureRDK(molecule):
	"""
	Given an RDKit molecule, create an OpenMM System and use to
	generate a ParmEd structure using the SMIRNOFF forcefield parameters.
	"""
	ff = pkg_resources.resource_filename('duck', "parameters/smirnoff99Frosst.offxml")
	with open(ff) as ffxml:
	    mol_ff = forcefield_rdk.ForceField(ffxml)
	#TODO : integrate charges
	# Can run this - should be already done
	# os.system("antechamber -i " + mol2_file + " -fi mol2 -o " + out_file + " -fo mol2 -at sybyl  -c bcc")
	charged_molecule = molecule
	mol_top, mol_sys, mol_pos = create_system_from_molecule_rdk(mol_ff, charged_molecule)
	molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
	return molecule_structure
