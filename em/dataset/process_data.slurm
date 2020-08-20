#!/bin/bash
#SBATCH --job-name=dataGen
#SBATCH --output=%j.%N.out
#SBATCH --nodes=12
#SBATCH --tasks-per-node=32
#SBATCH --partition=nu-wide
#SBATCH --ntasks=384	
#SBATCH --time=1-00:00:00

module load gcc/9.3.0
module load mpich/3.3.2-gcc-9.3.0 

source ~/.bashrc 
conda activate visual

time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 27 48
