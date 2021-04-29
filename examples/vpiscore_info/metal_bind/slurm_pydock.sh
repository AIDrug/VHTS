#!/bin/bash
#SBATCH --job-name=docking # Job name
#SBATCH --partition=docking # Job name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=256
#SBATCH --time=72:00:00 # Time limit hrs:min:sec
#SBATCH --output=subdock_%j.log # Standard output and error log

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
pydock_run.py --arg_file pydock_config.txt --my_module ./my_module.py

date
