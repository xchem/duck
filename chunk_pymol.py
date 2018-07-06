from pymol import cmd, stored
import os

ligand = "ligand"
protein = "protload"
mol_file = "MURD-x0349.mol"
prot_file = "MURD-x0349_apo.pdb"
out_save = "protein_out.pdb"
prot_protein_chunk = "protein_chunk_prot.pdb"
out_obj = "out_sele"
cutoff = 8

def return_tleap(out_save,prot_protein_chunk):
    return """source /opt/conda/envs/prepare/dat/leap/cmd/leaprc.ff14SB.redq
mol = loadpdb """+out_save+"""
savepdb mol """+prot_protein_chunk+"""
quit"""

# Load and select chunk 
cmd.load(mol_file,ligand)
cmd.load(prot_file,protein)
cmd.select(out_obj,"br. " + protein + "  within " + str(cutoff) + " of " + ligand)
cmd.select("no_het",out_obj + " and polymer")
cmd.save(out_save,"no_het")
# Now prepare with TLEAP
# Now do tleap
out_f = open("run.tleap","w")
out_f.write(return_tleap(out_save,prot_protein_chunk))
out_f.close()
os.system("/opt/conda/envs/prepare/bin/tleap -f run.tleap")
