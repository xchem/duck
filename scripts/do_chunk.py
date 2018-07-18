from duck.steps.chunk import chunk_with_amber,do_tleap,remove_prot_buffers_alt_locs,find_disulfides
import sys

def main(prot_file, mol_file, cutoff):
    # A couple of file name
    orig_file = prot_file
    chunk_protein = "protein_out.pdb"
    chunk_protein_prot = "protein_out_prot.pdb"
    # Do the removal of buffer mols and alt locs
    prot_file = remove_prot_buffers_alt_locs(prot_file)
    # Do the chunking and the protonation
    # Chunk
    chunk_with_amber(mol_file, prot_file, chunk_protein, cutoff)
    # Protontate
    disulfides = find_disulfides(chunk_protein)
    do_tleap(chunk_protein, chunk_protein_prot, disulfides)


if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],int(sys.argv[3]))