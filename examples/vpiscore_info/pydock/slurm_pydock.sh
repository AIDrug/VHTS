#!/bin/bash
#SBATCH --job-name=docking # Job name
#SBATCH --partition=docking # Job name
#SBATCH --ntasks=256 # Run on a single CPU
#SBATCH --time=72:00:00 # Time limit hrs:min:sec
#SBATCH --output=subdock_%j.log # Standard output and error log
## --nodes=1
## --ntasks-per-node=256

## --ntasks=256 # Run on a single CPU
## --mail-user=ned@aitrics.com
## --mail-type=All

pwd; hostname; date

NPROCS=$SLURM_NPROCS
echo NPROCS=$NPROCS
#NPROCS=$SLURM_NTASKS
#echo NPROCS=$NPROCS
#NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
#echo NPROCS=$NPROCS

source activate ml
conda deactivate
source activate ml


pydock_run.py --arg_file pydock_config.txt --my_module ./my_modlue.py --num_sub_proc $NPROCS
date
