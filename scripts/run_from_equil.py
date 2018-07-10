from duck.utils.get_from_fragalysis import get_from_prot_code
from duck.utils.cal_ints import find_interaction
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md
# Define these
prot_code = "MURD-x0374"
prot_int = "A_LYS_311_N"
# Cutoff for the sphere
# Now it does the magic
results  = get_from_prot_code(prot_code)
prot_file = results[0]
# Now find the interaction and save to a file

results = find_interaction(prot_int,prot_file)
startdist = results[2]
for i in range(20):
    if i==0:
        md_start = "equil.chk"
    else:
        md_start = "md_"+str(i-1)+".chk"
    perform_md(md_start,"md_"+str(i)+".chk","md_"+str(i)+".csv","md_"+str(i)+".pdb",md_len=0.5)
    # Now find the interaction and save to a file
    run_steered_md(300,"md_"+str(i)+".chk","smd_"+str(i)+"_300.csv","smd_"+str(i)+"_300.dat",
                   "smd_"+str(i)+"_300.pdb","smd_"+str(i)+"_300.dcd",startdist)
    run_steered_md(325,"md_"+str(i)+".chk","smd_"+str(i)+"_325.csv","smd_"+str(i)+"_325.dat",
                   "smd_"+str(i)+"_325.pdb","smd_"+str(i)+"_325.dcd",startdist)
