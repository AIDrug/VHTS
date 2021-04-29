
export PYTHONPATH="${YOUR_PATH}/VHTS:$PYTHONPATH"
export PATH="${YOUR_PATH}/VHTS/bin:$PATH"
# PIFinder
export PYTHONPATH="${YOUR_PATH}/PIFinder:$PYTHONPATH
export PATH="{YOUR_PATH}/PIFinder/bin:$PATH"

# prepare 
find_pf_receptor.py -v config.txt -r pdb/2G1MA_receptor.pdb -p pdb/pf_receptor.txt -t pdb/2G1MA_4HG.pdb -u pdb/pf_receptor_info.txt

# fix pdb/pf_receptor_info.txt to below (only remain metal and change weight)
vi pdb/pf_receptor_info.txt

$ echo pdb/pf_receptor_info.txt
feature_type:atom_idx_group:pattern_idx:pseudo_atom_coor:atom_type:etc:chain_id:residue_name:residue_num:weight
Metal:(3315):3:(39.864,22.370,13.612):Fe::A:FE2:404:5.000


### VHTS for batch system (slurm)

mkdir workflow
cd workflow
mkdir current done master todo
cd ..

cp smiles_list.csv workflow/master/remain.csv

# master
nohup master_dock.py --arg_file master_config.txt > master_log.txt 2>&1 &
# sub_dock 
sbatch slurm_sub_dock.sh
# 


### non VHTS 
pydock_run.py --arg_file pydock_config.txt --my_module ./my_module.py

