import os,glob,yaml
run_files = glob.glob("*/*/run.yaml")

checkpoints = {"ligand":"ANTECHAMBER_AM1BCC.AC","solvated": "complex_solvated.pdb", "merged": "complex.pdb",
     "equil": "equil.chk", "heating": "heating.csv","density": "density.csv","md_$i": "md_$i.chk"}
check_list = ["ligand","solvated", "merged", "equil", "heating","density","md_$i"]

out_d = {}
for yaml in run_files:
    base_dir = os.path.dirname(yaml)
    run_name = os.path.split(base_dir)[-1]
    run_dict = yaml.parse(open(yaml).read())
    checkpoints_list = []
    for check in check_list:
        file_paths= []
        if "$i" in checkpoints.check:
            for i in range(int(run_dict["num_smd_cycles"])):
                check_p = checkpoints[check].replace("$i", str(i))
                file_p = os.path.join(base_dir,check_p)
                checkpoints_list.append(check_p)
                checkpoints_list.append((file_p,check_p))
        else:
            check_p = checkpoints[check]
            file_p = os.path.join(base_dir,check_p)
            checkpoints_list.append((file_p, check_p))
    results = []
    for check_p,file_p in checkpoints_list:
        if os.path.isfile(file_p):
            results.append((check_p,True))
        else:
            results.append((check_p,False))
    out_d[run_name] = results
print(out_d)