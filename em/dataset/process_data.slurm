#!/bin/bash
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=12
##SBATCH --tasks-per-node=32
##SBATCH --partition=nu-wide
##SBATCH --ntasks=384	
##SBATCH --time=1-00:00:00

#SBATCH --job-name=proteins
#SBATCH --output=%j.%N.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=72
#SBATCH --partition=dribe
#SBATCH --ntasks=72   
#SBATCH --time=3-00:00:00


module load gcc/9.3.0
module load mpich/3.3.2-gcc-9.3.0 

source ~/.bashrc 
conda activate visual

#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 27 48
#probando generacion sintetica
time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 27 34 --g 30 --f 1
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 0 --f 1
