import sys
from vmd import evaltcl
import openbabel

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

# Define the input variables
prot_pdb = sys.argv[1]
prot_mol2 = prot_pdb.replace(".pdb",".mol2")
file_out_name=prot_pdb.replace(".pdb","out.mol2")
index = 1561#sys.argv[2]
cutoff = 6 #sys.argv[3]

# Convert the PDB to mol2 
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb","mol2")
mol = openbabel.OBMol()
obConversion.ReadFile(mol, prot_pdb)
mol.AddHydrogens()
obConversion.WriteFile(mol, prot_mol2)

# Convert into a VMD command from within Python

out_f = open("run.tcl","w")
out_f.write(return_tcl(prot_mol2,file_out_name,str(cutoff),str(index)))
out_f.close()
evaltcl('play run.tcl')

