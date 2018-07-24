import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import parmed
from duck.steps.chunk import remove_prot_buffers_alt_locs

def clean_up_mol(mol_block):
    mol = Chem.MolFromMolBlock(mol_block,removeHs=False)
    return Chem.MolToMolBlock(AllChem.AddHs(mol,addCoords=True))

def clean_up_prot(pdb_file):
    output_file = remove_prot_buffers_alt_locs(pdb_file)



def get_from_prot_code(prot_code):
    out_prot = prot_code+"_apo.pdb"
    out_mol = prot_code+".mol"
    r = requests.get("https://fragalysis.apps.xchem.diamond.ac.uk/api/proteins/",params={"code":prot_code})
    result = r.json()["results"][0]
    pdb_url = result["pdb_info"]
    pdb_id = result["id"]
    r = requests.get("https://fragalysis.apps.xchem.diamond.ac.uk/api/molecules/",params={"prot_id":pdb_id})
    out_json = r.json()["results"][0]
    sdf_info = clean_up_mol(out_json["sdf_info"])
    out_f = open(out_mol,"w")
    out_f.write(sdf_info)
    out_f.close()
    # Now s
    r = requests.get(pdb_url)
    out_f = open(out_prot,"w")
    out_f.write(r.text)
    out_f.close()
    return [out_prot,out_mol]

if __name__ == "__main__":
    prot_code = "MURD-x0373"
    get_from_prot_code(prot_code)