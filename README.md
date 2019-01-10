# Installation

Make a fresh Conda environment
```
git clone https://github.com/xchem/duck
cd duck
conda env create -f environment.yaml 
```
# Running

Activate conda and run like this:
```
source activate duck
frag_duck run.yaml
```

where run.yaml is a fiel like the following:

```
prot_code: '1n2v'
prot_int: 'A_ASP_156_OD2'
lig_id: 'BDI'
cutoff: 9
md_len: 0.5
distance: 2.5
init_velocity: 0.00001
num_smd_cycles: 1
gpu_id: '3'
apo_pdb_file: '1n2v_apo.pdb'
mol_file: ligand.mol
```
