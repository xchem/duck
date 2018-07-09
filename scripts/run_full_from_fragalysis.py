from duck.steps.parametrize import prepare_system
from duck.utils.get_from_fragalysis import get_from_prot_code
from duck.utils.cal_ints import find_interaction
from duck.steps.prepare import prep_lig
from duck.steps.chunk import chunk_with_pymol,do_tleap
from duck.steps.equlibrate import do_equlibrate
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md
# Define these
prot_code = "MURD-x0373"
prot_int = "A_LYS_311_N"
# Cutoff for the sphere
cutoff = 12
# A couple of file names
chunk_protein = "protein_out.pdb"
chunk_prot_protein = "protein_out_prot.pdb"

# Now it does the magic
results  = get_from_prot_code(prot_code)
prot_file = results[0]
mol_file = results[1]
# Chunk
chunk_with_pymol(mol_file, prot_file, chunk_protein, cutoff)
# Protontate
do_tleap(chunk_prot_protein,chunk_protein)
results = prep_lig(mol_file,prot_code)
mol2_file = results[0]
results = prepare_system(mol2_file,chunk_prot_protein)
complex = results[0]
# Now find the interaction and save to a file
results = find_interaction(prot_int,prot_file)
startdist = results[2]
# Now do the equlibration
results = do_equlibrate()
perform_md(results[0],"md_1.chk","md_1.csv","md_1.pdb")
# Now find the interaction and save to a file
results = find_interaction(prot_int,prot_file)
run_steered_md(300,"md_1.chk","smd_1_300.csv","smd_1_300.dat","smd_1_300.pdb","smd_1_300.dcd",startdist)
run_steered_md(325,"md_1.chk","smd_1_325.csv","smd_1_325.dat","smd_1_325.pdb","smd_1_325.dcd",startdist)

