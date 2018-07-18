from duck.steps.chunk import chunk_with_amber,do_tleap,remove_prot_buffers_alt_locs,find_disulfides
import sys,yaml

def main(prot_file, mol_file, interaction, cutoff):
    # A couple of file name
    orig_file = prot_file
    chunk_protein = "protein_out.pdb"
    chunk_protein_prot = "protein_out_prot.pdb"
    # Do the removal of buffer mols and alt locs
    prot_file = remove_prot_buffers_alt_locs(prot_file)
    # Do the chunking and the protonation
    # Chunk
    chunk_with_amber(mol_file,prot_file,interaction,chunk_protein,cutoff,orig_file)
    # Protontate
    disulfides = find_disulfides(chunk_protein)
    do_tleap(chunk_protein, chunk_protein_prot, disulfides)

if __name__ == "__main__":
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
    prot_file = out_data["apo_pdb_file"]
    mol_file = out_data["mol_file"]
    main(prot_file,mol_file,prot_int,cutoff)