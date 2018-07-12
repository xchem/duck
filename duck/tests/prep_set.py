import requests,sys
from rdkit import Chem
import parmed
from parmed.geometry import distance2
import operator
import math,os

def gen_yaml(pdb_code, interaction, lig_id, pdb_file, mol_file,gpu_id):
    return """prot_code: '"""+pdb_code+"""'
prot_int: '"""+interaction+"""'
lig_id: '"""+lig_id+"""'
cutoff: 12
md_len: 0.5
distance: 2.5
init_velocity: 0.00001
num_smd_cycles: 20
gpu_id: '"""+str(gpu_id)+"""'
apo_pdb_file: '"""+pdb_file+"""'
mol_file: '"""+mol_file+"""'
"""

def proc_dir(dir_name="1b9v_GLU276_OE2_3100",batch=0,gpu_id=0,new_batch=False):
    # This needs fixing for sure
    pdb_code = dir_name.split("_")[0].split("-")[0]
    res_name = dir_name.split("_")[1][:3]
    res_num = int(dir_name.split("_")[1][3:])
    atom_name = dir_name.split("_")[2]
    atom_num = dir_name.split("_")[3]
    chain_id = "A"
    batch_dir = "BATCH_"+str(batch)
    if new_batch == True:
        if batch != 0:
            os.chdir("../")
        if not os.path.isdir(batch_dir):
            os.mkdir(batch_dir)
        os.chdir(batch_dir)
    if os.path.isdir(dir_name):
        return
    os.mkdir(dir_name)
    os.chdir(dir_name)
    pdb_files = get_pdb_apo(pdb_code)
    lig_id = find_lig_id(pdb_files[0],pdb_files[2],res_name, res_num, atom_name)
    mol_file = get_ligand(lig_id,pdb_code)
    interaction = "_".join([str(x) for x in [chain_id,res_name,res_num,atom_name]])
    out_yaml = gen_yaml(pdb_code, interaction ,lig_id, pdb_files[0], mol_file,gpu_id)
    out_f = open("run.yaml","w")
    out_f.write(out_yaml)
    os.chdir("../")


def find_lig_id(pdb_file, het_pdb, res_name, res_num, atom_name):
    print(pdb_file)
    protein = parmed.load_file(pdb_file)
    atoms = []
    for residue in protein.residues:
        if residue.name == res_name:
            if residue.number == res_num:
                for atom in residue.atoms:
                    if atom.name == atom_name:
                        atoms.append(atom)
    if len(atoms)==0:
        print("ISSUE FINDING ATOM: "+" ".join([pdb_file,res_name, str(res_num), atom_name]))
        sys.exit()
    else:
        prot_atom = atoms[0]
    ligands = parmed.load_file(het_pdb)["!(:HOH,NA,CL)"]
    distance_atom_2 = [(x.residue.name, math.sqrt(distance2(x, prot_atom)),x.name) for x in ligands.atoms]
    distance_atom_2.sort(key=operator.itemgetter(1))
    print(distance_atom_2[0],distance_atom_2[1],distance_atom_2[2])
    return distance_atom_2[0][0]


def get_ligand(lig_id,pdb_code):
    url = "https://www.rcsb.org/pdb/download/downloadLigandFiles.do?ligandIdList="+lig_id+"&structIdList="+pdb_code+"&instanceType=all&excludeUnobserved=false&includeHydrogens=false"
    # Save
    mol_file = pdb_code + "_" +lig_id + ".mol"
    request = requests.get(url)
    mol = Chem.MolFromMolBlock(request.text)
    Chem.MolToMolFile(mol, mol_file)
    return mol_file

def write_pdb(lines,file_name):
    pdb_file = open(file_name,"w")
    for line in lines:
        pdb_file.write(line + "\n")
    pdb_file.close()

def get_pdb_apo(pdb_code):
    # Get the APO PDB (remove all HETATM records)
    url = "https://files.rcsb.org/download/"+pdb_code+".pdb"
    request = requests.get(url)
    all_lines = request.text.split("\n")
    apo_lines = [x for x in all_lines if x.startswith("ATOM")]
    het_lines = [x for x in all_lines if x.startswith("HETATM")]
    # Write the pdb files otu
    apo_pdb = pdb_code + "_apo.pdb"
    holo_pdb = pdb_code + ".pdb"
    het_pdb = pdb_code + "_het.pdb"
    write_pdb(apo_lines,apo_pdb)
    write_pdb(het_lines,het_pdb)
    write_pdb(all_lines,holo_pdb)
    return [apo_pdb,holo_pdb,het_pdb]


if __name__ == "__main__":
    dirs = ['2mcp_ARG52_NH1_781', '1yv3_LEU262_O_3916', '1q1g_GLU184_OE2_2791', '1n2v_ASP156_OD2_2255', '1lpz_ASP189_OD1_3562', '1l7f_GLU227_OE1_2279', '1k3u_SER235_OG_3501', '1hgi_GLY135_O_1984', '1hgh-1_GLY135_O_1985', '1fh8_GLN203_OE1_3004', '1f0t_ASP189_OD1_2470', '1exa_SER289_OG_1722', '1c1b_LYS101_O_1593', '1fh8_TRP273_NE1_4102', '1ydr_VAL123_N_1823', '1v0p_ASP85_OD1_1410', '1ulb_GLU201_OE1_3136', '1rob_PHE120_N_1789', '1n46_ARG316_NH1_1636', '1l7f_ARG118_NH1_568', '1l2s_ALA318_O_4782', '1k1j_ASP189_OD2_2471', '1ivd_ARG118_NH1_569', '1hwi_LYS691_NZ_9397', '1hgj_ASN137_N_1999', '1hgi_TYR98_OH_1438', '1hgh-1_SER136_OG_1992', '1fhd-2_GLU127_OE2_1871', '1fh9_GLU233_OE2_3481', '1fh9_ASN126_ND2_1856', '1f0u_SER190_OG_2481', '1exa_MET272_SD_1437', '1br6_VAL81_O_1295', '1b9v_GLU276_OE2_3100']
    old_batch = -1
    for i,dir in enumerate(dirs):
        print(dir)
        gpu_id = i % 4
        batch = i // 4
        print(batch,gpu_id)
        new_batch = batch!=old_batch
        print(new_batch)
        proc_dir(dir,batch,gpu_id,new_batch)
        old_batch=batch