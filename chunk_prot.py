import sys
from vmd import evaltcl
import openbabel
import re
import string
import os


def rename_residues(file_mol2,new_output_file):
    #Extrction of key from *.mol2 file for further translation
    key = {}
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
    g = open(new_output_file,'w')
    for line in f:
        a = line.strip()
        if a[0:4] == 'ATOM':
            index = int(re.sub("\D", " ", a[23:26]))
            g.write(a[:23]+res_chunk_pdb_mod[index-1]+a[26:]+'\n')
        else:
            g.write(a+'\n')
    f.close()
    g.close()


def return_rdock_tcl(prot_mol2, tmp2_file_out_name, temp_txt, index, cutoff):
    return '''###Arguments given from a command line are:
###vmd -dispdev text -e makechunk.tlc -args in_file_name out_file_name cavity_file index_ref_at cutoff

set file_name "'''+prot_mol2+'''"
set file_out "'''+tmp2_file_out_name+'''"
set cavity_file '''+temp_txt+'''

set ref_at_0 '''+index+'''
#renumerate because vmd starts with 0
set ref_at [expr $ref_at_0 - 1]
set cutoff '''+cutoff+'''

mol new $file_name


#selecting core chunk of residues that exist within cavity
set fp [open "$cavity_file" r]
set file_data [read $fp]
close $fp
set temp [split $file_data "\n"]
set cavity_id_0 [lreplace $temp end end]
set chunk_0_1 [atomselect 0 "resid $cavity_id_0"]

set chunk_0_1_ind [$chunk_0_1 get resid]

# adding residues within distace of 6A from main atom

set chunk_0_2 [atomselect 0 "protein within $cutoff of index $ref_at"]

set chunk_0_2_ind [$chunk_0_2 get resid]


#joining two selections

set chunk [atomselect 0 "resid $chunk_0_1_ind $chunk_0_2_ind"]


#caps is a selection of indexes of caps
set caps [atomselect 0 "none"]

#extractling list of residue identificantion numbers

set cavity_id_0 [$chunk get resid]
set cavity_id [lsort -unique -integer $cavity_id_0]

set protein_whole_0 [atomselect 0 "protein"]
set protein_whole_0_0 [$protein_whole_0 get resid]
set protein_whole [lsort -unique -integer $protein_whole_0_0]
set min [lindex $protein_whole 0]
set max [lindex $protein_whole end]

#selecting atoms from  +/-1 residue in order to create protecting groups

foreach i $cavity_id {


        set res_tmp [atomselect 0 "name CA and resid $i"]
        set res_num_tmp [$res_tmp get resname]
        regsub {^[A-Z]+} $res_num_tmp "" res_num_tmp2_0
        regsub {[A-Z]} $res_num_tmp2_0 "" res_num_tmp2
        set res_num_tmp2_min1 [expr $res_num_tmp2 - 1]
        set res_num_tmp2_plu1 [expr $res_num_tmp2 + 1]

        set res_minus_1 [expr $i - 1]

        set res_min_1_tmp [atomselect 0 "name CA and resid $res_minus_1"]
        set res_num_min_1_tmp [$res_min_1_tmp get resname]
        regsub {^[A-Z]+} $res_num_min_1_tmp "" res_num_min_1_tmp2_0
        regsub {[A-Z]} $res_num_min_1_tmp2_0 "" res_num_min_1_tmp2


        set res_plus_1 [expr $i + 1]

        set res_plu_1_tmp [atomselect 0 "name CA and resid $res_plus_1"]
        set res_num_plu_1_tmp [$res_plu_1_tmp get resname]
        regsub {^[A-Z]+} $res_num_plu_1_tmp "" res_num_plu_1_tmp2_0
        regsub {[A-Z]} $res_num_plu_1_tmp2_0 "" res_num_plu_1_tmp2



	if { $res_num_min_1_tmp2 != $res_num_tmp2_min1 || $i == $min } { set res_minus_1 [expr $i - 0] } \
	else { set res_minus_1 [expr $i - 1] }
	set cap_minus_1 [atomselect 0 "name CA C O and resid $res_minus_1"]
	set cap_minus_ind [$cap_minus_1 get index]

	if { $res_num_plu_1_tmp2 != $res_num_tmp2_plu1 || $i == $max } { set res_minus_1 [expr $i + 0] } \
	else { set res_plus_1 [expr $i + 1] }
	set cap_plus_1 [atomselect 0 "name CA N and resid $res_plus_1"]
	set cap_plus_ind [$cap_plus_1 get index]

	set caps_ind [$caps get index]
	set caps_temp [atomselect 0 "index $cap_plus_ind $cap_minus_ind $caps_ind"]
	set caps_temp_ind [$caps_temp get index]
	set caps [atomselect 0 "index $caps_temp_ind"]

}

#Selectin additional residues if gaps are less than 3 residues

set gaps [atomselect 0 "none"]

set var [lindex $cavity_id 0]

foreach i $cavity_id {

	if {$var == [expr $i - 2]} {
	set val2 [expr $i - 1]
	set gaps_temp [atomselect 0 "resid $val2"]
        set gaps_temp_id [$gaps_temp get resid]
        set gaps_id [$gaps get resid]
        set gaps [atomselect 0 "resid $gaps_temp_id $gaps_id"]

	} elseif {$var == [expr $i - 3]} {
	set val2 [expr $i - 1]
	set val3 [expr $i - 2]
	set gaps_temp [atomselect 0 "resid $val2 $val3"]
        set gaps_temp_id [$gaps_temp get resid]
	set gaps_id [$gaps get resid]
        set gaps [atomselect 0 "resid $gaps_temp_id $gaps_id "]

	}

	set var $i
}

#joining three selections
set index1 [$chunk get index]
set index2 [$caps get index]
set index3 [$gaps get index]
set chunk_final [atomselect 0 "index $index1 $index2 $index3 and not hydrogen"]
#set chunk_final [atomselect 0 "index $index1 $index2 and not hydrogen"]

#saving selected residues
$chunk_final writepdb $file_out
quit
'''


def return_tcl(prot_in,prot_out,cut_off,atom_index):
  return '''###Arguments given from a command line are:
###vmd -dispdev text -e make_chunk_2.tlc -args in_file_name_mol2 out_file_name.pdb index_ref_at cutoff

set file_name "'''+prot_in+'''"
set file_out "'''+prot_out+'''"
set ref_at_0 '''+atom_index+'''
#renumerate because vmd starts with 0
set ref_at [expr $ref_at_0 - 1]
set cutoff '''+cut_off+'''
mol new $file_name

# adding residues within distace of 6A from main atom

set chunk_0_2 [atomselect 0 "protein within $cutoff of index $ref_at"]

set chunk_0_2_ind [$chunk_0_2 get resid]

#joining two selections

set chunk [atomselect 0 "resid $chunk_0_2_ind"]

#caps is a selection of indexes of caps
set caps [atomselect 0 "none"]

#extractling list of residue identificantion numbers

set cavity_id_0 [$chunk get resid]
set cavity_id [lsort -unique -integer $cavity_id_0]


set protein_whole_0 [atomselect 0 "protein"]
set protein_whole_0_0 [$protein_whole_0 get resid]
set protein_whole [lsort -unique -integer $protein_whole_0_0]
set min [lindex $protein_whole 0]
set max [lindex $protein_whole end]


#selecting atoms from  +/-1 residue in order to create protecting groups


foreach i $cavity_id {

	set res_tmp [atomselect 0 "resid $i and name CA"]
	set res_num_tmp [$res_tmp get resname]


	regsub {^[A-Z]+} $res_num_tmp "" res_num_tmp2_0
        regsub {[A-Z]} $res_num_tmp2_0 "" res_num_tmp2
	set res_num_tmp2_min1 [expr $res_num_tmp2 - 1]
	set res_num_tmp2_plu1 [expr $res_num_tmp2 + 1]

	set res_minus_1 [expr $i - 1]
	set res_min_1_tmp [atomselect 0 "name CA and resid $res_minus_1"]
        set res_num_min_1_tmp [$res_min_1_tmp get resname]
        regsub {^[A-Z]+} $res_num_min_1_tmp "" res_num_min_1_tmp2_0
        regsub {[A-Z]} $res_num_min_1_tmp2_0 "" res_num_min_1_tmp2

        set res_plus_1 [expr $i + 1]
	set res_plu_1_tmp [atomselect 0 "name CA and resid $res_plus_1"]
        set res_num_plu_1_tmp [$res_plu_1_tmp get resname]
        regsub {^[A-Z]+} $res_num_plu_1_tmp "" res_num_plu_1_tmp2_0
        regsub {[A-Z]} $res_num_plu_1_tmp2_0 "" res_num_plu_1_tmp2

	if { $res_num_min_1_tmp2 != $res_num_tmp2_min1 || $i == $max } { set res_minus_1 [expr $i - 0] } \
	else { set res_minus_1 [expr $i - 1] }
	set cap_minus_1 [atomselect 0 "name CA C O and resid $res_minus_1"]
	set cap_minus_ind [$cap_minus_1 get index]

	if { $res_num_plu_1_tmp2 != $res_num_tmp2_plu1 || $i == $max } { set res_minus_1 [expr $i + 0] } \
	else { set res_plus_1 [expr $i + 1] }
	set cap_plus_1 [atomselect 0 "name CA N and resid $res_plus_1"]
	set cap_plus_ind [$cap_plus_1 get index]

	set caps_ind [$caps get index]
	set caps_temp [atomselect 0 "index $cap_plus_ind $cap_minus_ind $caps_ind"]
	set caps_temp_ind [$caps_temp get index]
	set caps [atomselect 0 "index $caps_temp_ind"]

}

#Selectin additional residues if gaps are less than 3 residues

set gaps [atomselect 0 "none"]

set var [lindex $cavity_id 0]

foreach i $cavity_id {

	if {$var == [expr $i - 2]} {
	set val2 [expr $i - 1]
	set gaps_temp [atomselect 0 "resid $val2"]
        set gaps_temp_id [$gaps_temp get resid]
        set gaps_id [$gaps get resid]
        set gaps [atomselect 0 "resid $gaps_temp_id $gaps_id"]

	} elseif {$var == [expr $i - 3]} {
	set val2 [expr $i - 1]
	set val3 [expr $i - 2]
	set gaps_temp [atomselect 0 "resid $val2 $val3"]
        set gaps_temp_id [$gaps_temp get resid]
	set gaps_id [$gaps get resid]
        set gaps [atomselect 0 "resid $gaps_temp_id $gaps_id "]

	}

	set var $i
}

#joining three selections
set index1 [$chunk get index]
set index2 [$caps get index]
set index3 [$gaps get index]
#set chunk_final [atomselect 0 "index $index1 and not hydrogen"]
#set chunk_final [atomselect 0 "index $index1 $index2 and not hydrogen"]
set chunk_final [atomselect 0 "index $index1 $index2 $index3 and not hydrogen"]

#saving selected residues
$chunk_final writepdb $file_out
quit    '''



import numpy as np
import math
import sys

# Gram-Schmidt process
def proj(v1, v2):  # projection fo vector v2 on vector v1
    return (np.dot(v2, v1) / np.dot(v1, v1)) * v1


def gs(V, X):  # return vector V moddified by orthogonalization to pair of vectors in X list
    u0 = X[0]
    u1 = X[1] - proj(X[0], X[1])
    return V - proj(u0, V) - proj(u1, V)


# function that returns triangle surface from atom list
def triangles(atoms):
    triangles = []
    n_at0 = 1
    for at0 in atoms:

        n_at1 = n_at0 + 1

        for at1 in atoms[n_at0:]:

            ref1 = at1 - at0

            for at2 in atoms[n_at1:]:

                ref2 = at2 - at0
                ref = [ref1, ref2]
                n_cross = 1

                for cross1 in atoms:

                    cross_ref1 = cross1 - at0
                    orthogonal_cross_ref1 = gs(cross_ref1, ref)

                    for cross2 in atoms[n_cross:]:

                        cross_ref2 = cross2 - at0
                        orthogonal_cross_ref2 = gs(cross_ref2, ref)

                        result = np.dot(orthogonal_cross_ref1, orthogonal_cross_ref2)

                        if result < -1.0e-8:
                            break

                    if result < -1.0e-8:
                        break

                    n_cross += 1

                if result >= -1.0e-8:
                    triangles.append([at0, at1, at2])

            n_at1 += 1
        n_at0 += 1
    return triangles


def get_atoms(file_pdb):  # gets atoms from pdb file
    atoms = []
    f = open(file_pdb, 'r')
    for line in f:
        if line[0:4] == 'ATOM':
            vect = np.array([float(line[31:39]), float(line[39:47]), float(line[47:55])])
            atoms.append(vect)

    f.close()
    return atoms


def get_bonds_mol2(file_mol2):  # get bond coordinates and assigns them to residue from mol2 file
    file_txt = []
    f = open(file_mol2)
    for line in f:
        file_txt.append(line.strip().split())

    f.close()
    try:
        start_atom = file_txt.index(['@<TRIPOS>ATOM']) + 1
        stop_atom = file_txt.index(['@<TRIPOS>BOND'])
        start_bond = stop_atom + 1
        stop_bond = file_txt.index(['@<TRIPOS>SUBSTRUCTURE'])
    except ValueError:
        return None

    bonds = []
    for line in file_txt[start_bond: stop_bond]:
        at1 = line[1]
        at2 = line[2]
        flag_hydrogen = 0
        for line_ref in file_txt[start_atom: stop_atom]:

            if line_ref[0] == at1:
                if line_ref[1][0] == 'H':
                    flag_hydrogen = 1
                    break
                coor_at1 = np.array([float(line_ref[2]), float(line_ref[3]), float(line_ref[4])])
                iden_at1 = line_ref[6]
            # else:
            #    continue
            if line_ref[0] == at2:
                if line_ref[1][0] == 'H':
                    flag_hydrogen = 1
                    break
                coor_at2 = np.array([float(line_ref[2]), float(line_ref[3]), float(line_ref[4])])
                iden_at2 = line_ref[6]
                break
        if flag_hydrogen == 0:
            bonds.append([coor_at1, coor_at2, iden_at1, iden_at2])

    return bonds



def narrow_sel2(atoms):  # narrows the selection to only surface atoms with vector method
    list_vect_tot = []

    for at1 in atoms:
        vect_tot = np.array([0.0, 0.0, 0.0])
        for at2 in atoms:
            vect_tot += (at2 - at1)

        dist_vect_tot = math.sqrt(np.dot(vect_tot, vect_tot))
        list_vect_tot.append([dist_vect_tot, at1])

    list_vect_tot.sort()

    list_vect_tot_short = list_vect_tot[int(len(list_vect_tot) * (1.0 - 0.5)):]  # i select

    list_vect_tot_short.reverse()

    spheres = []
    for at in list_vect_tot_short:
        spheres.append(at[1])

    atoms_2 = []
    for atom1 in spheres:
        flag = 1
        for atom2 in atoms_2:
            a = atom1 - atom2
            distance = math.sqrt(np.dot(a, a))
            if distance < 6:  # distance cutoff for group definition
                flag = 0
                break

        if flag == 1:
            atoms_2.append(atom1)
    return atoms_2


def intersec(triang, bond_list):
    intersections = []
    for triang_i in triang:
        V0 = triang_i[0]
        V1 = triang_i[1]
        V2 = triang_i[2]
        u = V1 - V0
        v = V2 - V0
        n = np.cross(u, v)
        n_norm = n / math.sqrt(np.dot(n, n))
        for bond_i in bond_list:
            P0 = bond_i[0]
            P1 = bond_i[1]
            if np.dot(n_norm, (P1 - P0)) != 0:
                r = np.dot(n_norm, (V0 - P0)) / np.dot(n_norm, (P1 - P0))

            if r > 0 and r < 1:
                PI = P0 + r * (P1 - P0)
                w = PI - V0
                s = (np.dot(u, v) * np.dot(w, v) - np.dot(v, v) * np.dot(w, u)) / (
                np.dot(u, v) * np.dot(u, v) - np.dot(u, u) * np.dot(v, v))
                t = (np.dot(u, v) * np.dot(w, u) - np.dot(u, u) * np.dot(w, v)) / (
                np.dot(u, v) * np.dot(u, v) - np.dot(u, u) * np.dot(v, v))

                if s > 0 and t > 0 and s + t < 1:
                    intersections.append(bond_i)
    return intersections

def surface_chunk_intesection(file_in_chunk_pdb,file_in_protein_mol2,file_out_list):
    ###Main body of the program

    atoms = get_atoms(file_in_chunk_pdb)
    bond_list = get_bonds_mol2(file_in_protein_mol2)
    if bond_list:
        atoms_2 = narrow_sel2(atoms)
        triang = triangles(atoms_2)
        intersections = intersec(triang, bond_list)
        list_add_res = []
        for i in intersections:
            if i[2] not in list_add_res:
                list_add_res.append(i[2])

            if i[3] not in list_add_res:
                list_add_res.append(i[3])
        output = open(file_out_list, 'w')
        for i in list_add_res:
            output.write(i + '\n')
        output.close()

        #########saving graphical outputs
        name_spheres = 'spheres_' + file_in_chunk_pdb.replace('.pdb', '.bild')
        name_triangle = 'triangle_' + file_in_chunk_pdb.replace('.pdb', '.bild')
        name_intersec = 'intersec_' + file_in_chunk_pdb.replace('.pdb', '.bild')
        name_bone = 'bone_' + file_in_chunk_pdb.replace('.pdb', '.bild')

        color = 'green'

        output = open(name_spheres, 'w')
        for i in atoms_2:
            output.write('.color ' + color + '\n')
            output.write('.sphere ')
            for ii in i:
                output.write(str(ii) + ' ')
            output.write('1.0 \n')
        output.close()

        output = open(name_triangle, 'w')
        for i in triang:
            output.write('.color ' + color + '\n')
            output.write('.polygon ')
            for ii in i:
                for iii in ii:
                    output.write(str(iii) + ' ')
            output.write('\n')
        output.close()

        output = open(name_bone, 'w')
        for tri in triang:
            output.write('.color ' + color + '\n')
            output.write(
                '.cylinder ' + str(tri[0][0]) + ' ' + str(tri[0][1]) + ' ' + str(tri[0][2]) + ' ' + str(tri[1][0]) + ' ' + str(
                    tri[1][1]) + ' ' + str(tri[1][2]) + ' 0.1\n')
            output.write('.color ' + color + '\n')
            output.write(
                '.cylinder ' + str(tri[0][0]) + ' ' + str(tri[0][1]) + ' ' + str(tri[0][2]) + ' ' + str(tri[2][0]) + ' ' + str(
                    tri[2][1]) + ' ' + str(tri[2][2]) + ' 0.1\n')
            output.write('.color ' + color + '\n')
            output.write(
                '.cylinder ' + str(tri[2][0]) + ' ' + str(tri[2][1]) + ' ' + str(tri[2][2]) + ' ' + str(tri[1][0]) + ' ' + str(
                    tri[1][1]) + ' ' + str(tri[1][2]) + ' 0.1\n')
        output.close()

        color_2 = 'orange'
        output = open(name_intersec, 'w')
        for i in intersections:
            output.write('.color ' + color_2 + '\n')
            output.write('.cylinder ')
            for coor in i[:2]:
                for xyz in coor:
                    output.write(str(xyz) + ' ')
            output.write('0.25\n')
        output.close()
    else:
        output = open(file_out_list, 'w')
        output.close()



def split_chains_rename_caps(prot_one,out_file):
    # reading the file into array 2xN containing list of splited fields and raw line for modification

    matrix = []
    f = open(prot_one, 'r')
    for line in f:
        matrix.append(
            [line.strip(), line.strip()])  # firs one is for searching patterns and the second one for modification
    f.close()

    # extract the list of numbers of residues
    all_res = []
    residues = []
    for i in matrix:
        if i[0][0:4] == 'ATOM':
            all_res.append(i[0][22:26].translate(None, string.letters))
            if i[0][22:26].translate(None, string.letters) not in residues:
                residues.append(i[0][22:26].translate(None, string.letters))

    # extract list of ACE and NME residues:
    res_ACE = []
    res_NME = []
    res_ACE.append(residues[0])
    ref = res_ACE[0]

    for i in residues[1:]:
        if int(i) - 1 != int(ref):
            res_ACE.append(i)
            res_NME.append(ref)
        ref = i

    res_NME.append(ref)

    # listing residues that are not caps, to comfirm later. beacuse it might be that first/last residue in
    # chain is also terminal residue in the whole protein, so it shouldn't be named as cap

    Not_caps = []
    for i in residues:
        nr_i = all_res.count(i)
        if nr_i > 3:
            Not_caps.append(i)

    # changing files and saving them to output
    for i in matrix:
        if i[0][0:4] == 'ATOM':
            name_i = i[0][22:26].translate(None, string.letters)
            if name_i in res_ACE and name_i not in Not_caps:
                i[1] = i[1][:17] + 'ACE' + i[1][20:]
                if i[0][12:16] == ' CA ':
                    i[1] = i[1][:12] + ' CH3' + i[1][16:]


            elif name_i in res_NME and name_i not in Not_caps:
                i[1] = i[1][:17] + 'NME' + i[1][20:]
                if i[0][12:16] == ' CA ':
                    i[1] = i[1][:12] + ' CH3' + i[1][16:]

    # saving the lines and including TER at the enad of every chain
    f = open(out_file, 'w')
    for i in matrix:
        f.write(i[1] + '\n')
        if i[0][0:4] == 'ATOM':
            if i[0][22:26].translate(None, string.letters) in res_NME and i[1][12:16] == ' CH3':
                f.write('TER\n')

    f.close()


def return_tleap():
    return """source /opt/conda/envs/prepare/dat/leap/cmd/leaprc.ff14SB.redq
mol = loadpdb tmp.pdb
savepdb mol tmp_out.pdb
quit"""
# Define the input variables
prot_pdb = sys.argv[1]
prot_mol2 = prot_pdb.replace(".pdb",".mol2")
file_out_name=prot_pdb.replace(".pdb","out.mol2")
# TODO add code to turn coordinate into an index. Input will be two co-rdinates for the ligand and protein interactions
# TODO
index = 1561#sys.argv[2]
cutoff = 6 #sys.argv[3]

# Convert the PDB to mol2 
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb","mol2")
mol = openbabel.OBMol()
obConversion.ReadFile(mol, prot_pdb)
mol.AddHydrogens()
obConversion.WriteFile(mol, prot_mol2)

# Run the first TCL script
out_f = open("run.tcl","w")
out_f.write(return_tcl(prot_mol2,"tmp_"+file_out_name,str(cutoff),str(index)))
out_f.close()
evaltcl('play run.tcl')

#patching gaps in structure by selecting additional residues
surface_chunk_intesection("tmp_"+file_out_name,prot_mol2,"temp2.txt")

out_f = open("run_two.tcl","w")
out_f.write(return_rdock_tcl(prot_mol2,"tmp2_" + file_out_name, "temp2.txt", str(index), str(cutoff)))
out_f.close()
evaltcl('play run_two.tcl')

#spliting chains and renamich caps
split_chains_rename_caps("tmp2_" + file_out_name, "tmp.pdb")

# Now do tleap
out_f = open("run.tleap","w")
out_f.write(return_tleap())
out_f.close()
os.system("/opt/conda/envs/prepare/bin/tleap -f run.tleap")
# Rename the residues
rename_residues(prot_mol2, file_out_name)