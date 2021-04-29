find_pf_receptor.py -v config.txt -r pdb/receptor.pdb -p pdb/pf_receptor.txt -t pdb/crystal_ligand.pdb -u pdb/pf_receptor_info.txt

# vhts for slurm
mkdir workflow
cd workflow
mkdir current done master todo
cd ..
cp smiles_list.pkl workflow/master/remain.pkl
nohup master_dock.py --arg_file master_config.txt > mm.txt 2>&1 &
sbatch slurm_sub_dock.sh


