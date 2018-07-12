from simtk.openmm import app
import parmed
from duck.steps.chunk import chunk_with_amber,do_tleap,remove_prot_buffers_alt_locs


def test_prep(mol_file,prot_file):
    chunk_protein = "protein_out.pdb"
    chunk_protein_prot = "protein_out_prot.pdb"
    # Now it does the magic
    # Chunk
    removed = remove_prot_buffers_alt_locs(prot_file)
    chunk_with_amber(mol_file, removed, chunk_protein, 12)
    # Protontate
    do_tleap(chunk_protein, chunk_protein_prot)
    protein = parmed.load_file(chunk_protein_prot)
    protein.write_pdb("fixed.pdb")
    protein = parmed.load_file("fixed.pdb")["!(:HOH,NA,CL,SO4,EDO)"]
    forcefield_str = 'amber99sb.xml'
    forcefield = app.ForceField(forcefield_str)
    protein_system = forcefield.createSystem(protein.topology)
    protein_pmd = parmed.openmm.load_topology(protein.topology, protein_system, protein.positions)