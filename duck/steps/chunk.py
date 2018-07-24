import os,pkg_resources
from rdkit import Chem
import parmed
from parmed.geometry import distance2
from duck.utils.cal_ints import find_atom

def return_tleap(prot_protein_chunk, out_save, disulfides=[]):
    param_f_path = pkg_resources.resource_filename('duck', "parameters/tleap/leaprc.ff14SB.redq")
    bond_list = ["bond mol."+str(num[0])+".SG mol."+str(num[1])+".SG " for num in disulfides]
    return """source """+param_f_path+"""
mol = loadpdb """ + prot_protein_chunk + """
"""+"\n".join(bond_list)+"""
savepdb mol """ + out_save + """
quit"""


def do_tleap(prot_protein_chunk, out_save, disulfides=[]):
    # Now do tleap
    out_f = open("run.tleap","w")
    out_f.write(return_tleap(prot_protein_chunk, out_save,disulfides))
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
    single_joins = {}
    double_joins = {}
    atom_set = set()
    for resid in residues:
        residue_set = set()
        double_res_set = set()
        for atom in resid.atoms:
            for bond_1 in atom.bonds:
                if bond_1.measure() >3.0:
                    print("Skipping too long bond (" + str(bond_1.measure()) + ") :  ", bond_1.atom1, bond_1.atom2)
                    continue
                if bond_1.atom1==atom:
                    x = bond_1.atom2
                else:
                    x = bond_1.atom1
                residue_set.add(x.residue)
                for atom_two in x.residue.atoms:
                    for bond_2 in atom_two.bonds:
                        if bond_2.measure() > 3.0:
                            print("Skipping too long bond ("+str(bond_2.measure())+") :  ", bond_2.atom1, bond_2.atom2)
                            continue
                        if bond_2.atom1 == atom:
                            conn_two = bond_2.atom2
                        else:
                            conn_two = bond_2.atom1
                        double_res_set.add(conn_two.residue)
                atom_set = add_cap(x,atom_set)
        double_res_set.difference_update(residue_set)
        single_joins[resid] = residue_set
        double_joins[resid] = double_res_set
    return single_joins,double_joins,atom_set

def find_neighbours(residues):
    single_joins, double_joins, atom_set = find_neighbour_residues(residues)
    new_residues = set()
    for resid_one in single_joins:
        for resid_two in single_joins:
            if resid_one == resid_two: continue
            single_join = single_joins[resid_two].intersection(single_joins[resid_one])
            # If it's just a terminal residue
            if resid_one not in single_joins[resid_two]:
                double_join = double_joins[resid_two].intersection(single_joins[resid_one])
            else:
                double_join = set()
            new_set = single_join.union(double_join)
            for new_res in new_set:
                residues.add(new_res)
                new_residues.add(new_res)
    return new_residues,atom_set

def convert_to_ace_nme(subset):
    remove_res_ids=[]
    remove_atom_ids = []
    for residue in subset.residues:
        if len(residue)==3:
            if set([x.name for x in residue.atoms]) == set(["CA","C","O"]):
                residue.name="ACE"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
            # If it's a proline
            if set([x.name for x in residue.atoms]) == set(["CA","CD","N"]):
                residue.name = "NME"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
                    if atom.name == "CD":
                        remove_atom_ids.append(str(atom.idx+1))
        elif len(residue)==2:
            if set([x.name for x in residue.atoms]) == set(["CA","N"]):
                residue.name="NME"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
            elif set([x.name for x in residue.atoms]) == set(["CB","SG"]):
                remove_res_ids.append(str(residue.idx+1))
    if remove_atom_ids != []:
        subset = subset["!(@"+",".join(remove_atom_ids)+")"]
    if remove_res_ids != []:
        subset = subset["!(:"+",".join(remove_res_ids)+")"]
    return subset


def remove_prot_buffers_alt_locs(prot_file):
    output_file = "no_buffer_altlocs.pdb"
    solvents = ["NA","CL","SO4","EDO","FMT","P04","DMS","EPE"]
    # Remove hydrogens and solvents and buffers
    protein = parmed.load_file(prot_file)["!(:HOH,"+",".join(solvents)+")"]["!(:=@H=)"]
    protein.write_pdb(output_file, altlocs="first")
    return output_file


def find_disulfides(input_file, threshold=6.2):
    structure = parmed.load_file(input_file)
    sulfurs = [x for x in structure.atoms if x.residue.name == "CYS" and x.name == "SG"]
    disulfides = []
    for atom_one  in sulfurs:
        for atom_two in sulfurs:
            if atom_one.idx >= atom_two.idx: continue
            dist = distance2(atom_one,atom_two)
            if dist < threshold:
                atom_one.residue.name = "CYX"
                atom_two.residue.name = "CYX"
                disulfides.append((atom_one.residue.number,atom_two.residue.number))
    structure.write_pdb(input_file)
    return disulfides

def find_res_idx(protein, chain, res_name, res_num):
    for residue in protein.residues:
        if residue.chain == chain:
            if residue.name==res_name:
                if residue.number==res_num:
                    return residue.idx + 1

def chunk_with_amber(mol_file="MURD-x0349.mol", prot_file="MURD-x0349_apo.pdb", interaction="A_LYS_311_N", out_save="protein_out.pdb", cutoff=9.0, orig_prot="MURD-x0349_apo.pdb"):
    # Load up the topology
    mol = Chem.MolFromMolFile(mol_file)
    pdb_mol_file = mol_file.replace(".mol",".pdb")
    Chem.MolToPDBFile(mol,pdb_mol_file)
    protein = parmed.load_file(prot_file)
    # get these details
    atom_idx,prot_atom = find_atom(interaction, orig_prot, protein)
    mask = parmed.amber.AmberMask(protein, "@"+str(atom_idx[0][0]+1)+"<:"+str(cutoff))
    residues = set([protein.atoms[i].residue for i, x in enumerate(mask.Selection()) if x == 1])
    residues = set([x for x in residues if x.name != "UNL"])
    # # Find all the residues that are connected to two residues in this list of residues.
    new_residues,atom_set = find_neighbours(residues)
    # Collect the atoms
    atom_idx = []
    for res in residues:
        atom_idx.extend([x.idx for x in res.atoms if x.altloc in ["","A"]])
    atom_idx.extend(atom_set)
    subset = protein[atom_idx]
    subset.write_pdb(out_save.replace(".pdb","_no_ace_nme.pdb"))
    subset = convert_to_ace_nme(parmed.load_file(out_save.replace(".pdb","_no_ace_nme.pdb")))
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