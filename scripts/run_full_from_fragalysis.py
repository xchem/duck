from duck.steps.parametrize import prepare_system
from duck.utils.get_from_fragalysis import get_from_prot_code
from duck.utils.cal_ints import find_interaction
from duck.steps.prepare import prep_lig
from duck.steps.chunk import chunk_with_amber,do_tleap,remove_prot_buffers_alt_locs,find_disulfides
from duck.steps.equlibrate import do_equlibrate
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md
from duck.utils.check_system import check_if_equlibrated
import yaml, sys, os

def run_simulation(prot_file, mol_file, prot_code, prot_int, cutoff, init_velocity, num_smd_cycles, gpu_id, md_len,params):
    if not os.path.isfile("equil.chk"):
        # A couple of file name
        orig_file = prot_file
        chunk_protein = "protein_out.pdb"
        chunk_protein_prot = "protein_out_prot.pdb"
        # Do the removal of buffer mols and alt locs
        if params.get("remove_buffers",True):
            prot_file = remove_prot_buffers_alt_locs(prot_file)
        # Do the chunking and the protonation
        if params.get("prep_chunk",True):
            # Chunk
            chunk_with_amber(mol_file, prot_file, prot_int, chunk_protein, cutoff, orig_file)
            # Protontate
            disulfides = find_disulfides(chunk_protein)
            do_tleap(chunk_protein,chunk_protein_prot,disulfides)
        else:
            chunk_protein_prot = prot_file
        # Paramaterize the ligand
        mol2_file = params.get("mol2_file_prepped", None)
        if not mol2_file:
            results = prep_lig(mol_file,prot_code)
            mol2_file = results[0]
        results = prepare_system(mol2_file,chunk_protein_prot)
        complex = results[0]
        # Now find the interaction and save to a file
        results = find_interaction(prot_int,orig_file)
        startdist = params.get("distance", results[2])
        # Now do the equlibration
        results = do_equlibrate(gpu_id=gpu_id)
    else:
        # Now find the interaction and save to a file
        results = find_interaction(prot_int,prot_file)
        startdist = params.get("distance", results[2])
    # Open the file and check that the potential is stable and negative
    if not check_if_equlibrated("density.csv",1):
        print("SYSTEM NOT EQUILIBRATED")
        sys.exit()
    if params.get("just_equilib",False):
        print("SYSTEM FENISHED EQULIBRATING")
        return
    # Now do the MD
    for i in range(num_smd_cycles):
        if i==0:
            md_start = "equil.chk"
        else:
            md_start = "md_"+str(i-1)+".chk"
        log_file = "md_"+str(i)+".csv"
        perform_md(md_start,"md_"+str(i)+".chk",log_file,"md_"+str(i)+".pdb",md_len=md_len,gpu_id=gpu_id)
        # Open the file and check that the potential is stable and negative
        if not check_if_equlibrated(log_file,3):
            print("SYSTEM NOT EQUILIBRATED")
            sys.exit()
        # Now find the interaction and save to a file
        run_steered_md(300,"md_"+str(i)+".chk","smd_"+str(i)+"_300.csv","smd_"+str(i)+"_300.dat",
                       "smd_"+str(i)+"_300.pdb","smd_"+str(i)+"_300.dcd",startdist,init_velocity=init_velocity,gpu_id=gpu_id)
        run_steered_md(325,"md_"+str(i)+".chk","smd_"+str(i)+"_325.csv","smd_"+str(i)+"_325.dat",
                       "smd_"+str(i)+"_325.pdb","smd_"+str(i)+"_325.dcd",startdist,init_velocity=init_velocity,gpu_id=gpu_id)


def main():
    # Define these in a YAML
    out_data = yaml.load(open(sys.argv[1]))
    prot_code = out_data["prot_code"]
    prot_int = out_data["prot_int"]
    # Cutoff for the sphere
    cutoff = out_data["cutoff"]
    md_len = out_data["md_len"]
    init_velocity = out_data["init_velocity"]
    num_smd_cycles = out_data["num_smd_cycles"]
    gpu_id = str(out_data["gpu_id"])
    # Now get the data from Fragalysis
    if "apo_pdb_file" not in out_data or "mol_file" not in out_data:
        results = get_from_prot_code(prot_code)
        prot_file = results[0]
        mol_file = results[1]
    else:
        prot_file = out_data["apo_pdb_file"]
        mol_file = out_data["mol_file"]
    run_simulation(prot_file, mol_file, prot_code, prot_int, cutoff, init_velocity, num_smd_cycles, gpu_id, md_len, out_data)

if __name__ == "__main__":
    main()
