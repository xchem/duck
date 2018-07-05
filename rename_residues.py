#script fixes numeration in chunk and creates the same one as in PDB file
import sys
import re

#Extrction of key from *.mol2 file for further translation
key = {}
file_mol2 = sys.argv[1]
check = []
f = open(file_mol2, 'r')
for line in f:
    a = line.strip().split()
    if len(a)>8:
        key_i = [int(re.sub("\D", " ", a[6])), int(re.sub("\D", " ", a[7]))]
        if key_i[0] not in check:
            check.append(key_i[0])
            key[str(key_i[0])] = str(key_i[1])
f.close()

#extracting list of residues to moddify
res_chunk_mol2 = []
f = open('tmp.pdb','r')
for line in f:
    a = line
    if a[:4] == 'ATOM' and int(re.sub("\D", " ", a[23:26])) not in res_chunk_mol2:
        res_chunk_mol2.append(int(re.sub("\D", " ", a[23:26])))
f.close()

#translating that list to old pdb numeration and moddifying format
res_chunk_pdb = []
for nr in res_chunk_mol2:
    res_chunk_pdb.append(key[str(nr)])

res_chunk_pdb_mod = []
for nr in res_chunk_pdb:
    if len(nr) == 1:
        res_chunk_pdb_mod.append('  '+nr)
    elif len(nr) == 2:
        res_chunk_pdb_mod.append(' '+nr)
    else:
        res_chunk_pdb_mod.append(nr)

#moddifying old file
f = open('tmp_out.pdb','r')
g = open(str(sys.argv[2]),'w')
for line in f:
    a = line.strip()
    if a[0:4] == 'ATOM':
        index = int(re.sub("\D", " ", a[23:26]))
        g.write(a[:23]+res_chunk_pdb_mod[index-1]+a[26:]+'\n')
    else:
        g.write(a+'\n')
f.close()
g.close()
