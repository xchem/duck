from pymol import cmd, stored
import os,pkg_resources


def return_tleap(out_save, prot_protein_chunk):
    param_f_path = pkg_resources.resource_filename('duck', "parameters/tleap/leaprc.ff14SB.redq")
    return """source """+param_f_path+"""
mol = loadpdb """ + out_save + """
savepdb mol """ + prot_protein_chunk + """
quit"""


def do_tleap(out_save, prot_protein_chunk):
    # Now do tleap
    out_f = open("run.tleap","w")
    out_f.write(return_tleap(out_save, prot_protein_chunk))
    out_f.close()
    os.system("tleap -f run.tleap")


def chunk_with_pymol(mol_file="MURD-x0349.mol", prot_file="MURD-x0349_apo.pdb", out_save="protein_out.pdb", cutoff=12):
    """
    :param mol_file:
    :param prot_file:
    :param cutoff:
    :return:
    """
    ligand = "ligand"
    protein = "protload"
    out_obj = "out_sele"
    # Load and select chunk
    cmd.load(mol_file,ligand)
    cmd.load(prot_file,protein)
    cmd.select(out_obj,"br. " + protein + "  within " + str(cutoff) + " of " + ligand)
    cmd.select("no_het",out_obj + " and polymer")
    cmd.save(out_save,"no_het")
    # Now do the capping
    # # add planer, trivalent nitrogen onto C terminus
    # edit cys////C
    # attach N,3,3
    # # edit that nitrogen
    # edit (elem N and neighbor cys////C)
    # # attach a tetrahedral methyl (note random position)
    # attach C,4,4
    # # here's an example of adding a whole residue from the library
    # edit cys////N
    # editor.attach_amino_acid("pk1","ace")

if __name__ == "__main__":
    chunk_with_pymol()
    do_tleap()