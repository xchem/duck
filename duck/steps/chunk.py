import os,pkg_resources
from rdkit import Chem
import parmed

def return_tleap(prot_protein_chunk, out_save):
    param_f_path = pkg_resources.resource_filename('duck', "parameters/tleap/leaprc.ff14SB.redq")
    return """source """+param_f_path+"""
mol = loadpdb """ + prot_protein_chunk + """
savepdb mol """ + out_save + """
quit"""


def do_tleap(prot_protein_chunk, out_save):
    # Now do tleap
    out_f = open("run.tleap","w")
    out_f.write(return_tleap(prot_protein_chunk, out_save))
    out_f.close()
    os.system("tleap -f run.tleap")

def add_cap(x, atom_set):
    # Add the cap
    atom_set.add(x.idx)
    # Add the capping atoms
    for y in x.bond_partners:
        atom_set.add(y.idx)
    return atom_set

def find_neighbour_residues(residues):
    out_d = {}
    atom_set = set()
    for resid in residues:
        residue_set = set()
        for atom in resid.atoms:
            for x in atom.bond_partners:
                residue_set.add(x.residue)
                atom_set = add_cap(x,atom_set)
        out_d[resid] = residue_set
    return out_d,atom_set

def find_neighbours(residues):
    out_d,atom_set = find_neighbour_residues(residues)
    for resid_one in out_d:
        for resid_two in out_d:
            if resid_one == resid_two: continue
            join = out_d[resid_two].intersection(out_d[resid_one])
            for new_res in join:
                residues.add(new_res)
    return residues,atom_set

def convert_to_ace_nme(subset):
    for residue in subset.residues:
        if len(residue)==3:
            if set([x.name for x in residue.atoms]) == set(["CA","C","O"]):
                residue.name="ACE"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
        elif len(residue)==2:
            if set([x.name for x in residue.atoms]) == set(["CA","N"]):
                residue.name="NME"
                print("HERE")
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
    return subset


def chunk_with_amber(mol_file="MURD-x0349.mol", prot_file="MURD-x0349_apo.pdb", out_save="protein_out.pdb", cutoff=7.0):
    # Load up the topology
    protein = parmed.load_file(prot_file)["!(:HOH,NA,CL,SO4,EDO)"]
    protein.write_pdb("no_altlocs.pdb", altlocs="first")
    protein = parmed.load_file("no_altlocs.pdb")
    mol = Chem.MolFromMolFile(mol_file)
    pdb_mol_file = mol_file.replace(".mol",".pdb")
    Chem.MolToPDBFile(mol,pdb_mol_file)
    ligand = parmed.load_file(pdb_mol_file)
    merged = ligand + protein
    mask = parmed.amber.AmberMask(merged, ":1<:"+str(cutoff))
    residues = set([merged.atoms[i].residue for i, x in enumerate(mask.Selection()) if x == 1])
    residues = set([x for x in residues if x.name != "UNL"])
    # # Find all the residues that are connected to two residues in this list of residues.
    neighbour_residues,atom_set = find_neighbours(residues)
    # Collect the atoms
    atom_idx = []
    for res in residues:
        atom_idx.extend([x.idx for x in res.atoms if x.altloc in ["","A"]])
    atom_idx.extend(atom_set)
    subset = merged[atom_idx]
    subset = convert_to_ace_nme(subset)
    subset.write_pdb(out_save)
    add_ter_records(out_save,out_save)
    return [out_save]

def prot_with_pdb_fixer(chunk_protein, chunk_prot_protein):
    os.system("pdbfixer " + chunk_protein + " --replace-nonstandard --output=" + chunk_prot_protein)
    return [chunk_prot_protein]

def add_ter_records(input_file,output_file):
    lines = open(input_file).readlines()
    output_f = open(output_file,"w")
    for line in lines:
        output_f.write(line)
        if "NME" in line:
            output_f.write("TER\n")
    return [output_file]


if __name__ == "__main__":
    chunk_with_amber()
    do_tleap()