import os,pkg_resources
from rdkit import Chem
import parmed

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


def chunk_with_amber(mol_file="MURD-x0349.mol", prot_file="MURD-x0349_apo.pdb", out_save="protein_out.pdb", cutoff=7.0):
    # Load up the topology
    protein = parmed.load_file(prot_file)["!(:HOH,NA,CL)"]
    mol = Chem.MolFromMolFile(mol_file)
    pdb_mol_file = mol_file.replace(".mol",".pdb")
    Chem.MolToPDBFile(mol,pdb_mol_file)
    ligand = parmed.load_file(pdb_mol_file)
    merged = ligand + protein
    mask = parmed.amber.AmberMask(merged, ":1<:"+str(cutoff))
    residues = set([merged.atoms[i].residue for i, x in enumerate(mask.Selection()) if x == 1])
    # Find all the residues that are connected to two residues in this list of residues.
    out_d = {}
    for resid in residues:
        residue_set = set()
        for atom in resid.atoms:
            for new_res in set([x.residue for x in atom.bond_partners]):
                residue_set.add(new_res)
        out_d[resid] = residue_set
    for resid_one in out_d:
        for resid_two in out_d:
            if resid_one == resid_two: continue
            join = out_d[resid_two].intersection(out_d[resid_one])
            if len(join) > 0:
                for new_res in join:
                    residues.add(new_res)
    residues = [x for x in residues if x.name != "UNL"]
    # Collect the atoms
    atom_idx = []
    for res in residues:
        atom_idx.extend([x.idx for x in res.atoms])
    merged[atom_idx].write_pdb(out_save)
    return [out_save]


if __name__ == "__main__":
    chunk_with_amber()
    do_tleap()