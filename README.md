[![Build Status](https://travis-ci.org/abradle/duck.svg?branch=master)](https://travis-ci.org/abradle/duck)
[![stable](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.1.0-blue.svg?style=flat)](https://github.com/abradle/duck)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/abradle/duck/blob/master/LICENSE.txt)

# Installation

## Conda

Make a fresh Conda environment
```
git clone https://github.com/abradle/duck
cd duck
conda env create -f environment.yaml 
```

## Docker

Pull down the latest image
```
docker pull abradle/duck
```

#### Running

Activate conda and run like this:
```
source activate duck
frag_duck run.yaml
```

Or with Docker run like this:
```
docker run -it -v $PWD:/data abradle/duck /bin/bash -c "frag_duck /data/run.yaml"
```

where run.yaml is a file like the following:

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
