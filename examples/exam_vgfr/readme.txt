# setup path
export PYTHONPATH="${YOUR_PATH}/vhts:$PYTHONPATH"
export PATH="${YOUR_PATH}/vhts/bin:$PATH"

# find docking box from template ligand file
cal_box.py --arg_file box_config.txt
# single docking 
pydock.py --arg_file pydock_config.txt --dock_config config.txt



# for batch system
# prepare 
mkdir workflow
cd workflow
mkdir current done master todo
cd ..

cp smiles_list.txt workflow/master/remain.txt
# master
nohup master_dock.py --arg_file master_config.txt > mm.txt 2>&1 &
# sub_dock 
nohup sub_dock.py --arg_file subdock_config.txt --dock_config config.txt --out_log file > a.txt 2>&1 &
nohup sub_dock.py --arg_file subdock_config.txt --dock_config config.txt --out_log file > b.txt 2>&1 &


