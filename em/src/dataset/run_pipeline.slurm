#!/bin/bash

#nu
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=25
##SBATCH --tasks-per-node=32
##SBATCH --partition=nu-all
##SBATCH --ntasks=800	
##SBATCH --time=1-00:00:00

#Dribe
#SBATCH --job-name=protpipe
#SBATCH --output=%j.%N.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=70
#SBATCH --partition=dribe
#SBATCH --ntasks=70
#SBATCH --time=3-00:00:00

#ANDALAN
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=1
##SBATCH --tasks-per-node=48
##SBATCH --partition=andalan
##SBATCH --ntasks=48
##SBATCH --time=24:00:00

#NUKWA
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=1
##SBATCH --tasks-per-node=48
##SBATCH --partition=nukwa
##SBATCH --ntasks=48    
##SBATCH --time=1-00:00:00
##SBATCH --nodelist=nukwa-04.cnca
module load gcc/9.3.0
module load mpich/3.3.2-gcc-9.3.0 

source ~/.bashrc 
conda activate visual

export PYTHONPATH=$PYTHONPATH:/work/mzumbado/biostructure

#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 22 29 --g 30 --f 1 --p 1
time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures feature_extraction.py

#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 1  --v 1 --s 22 29 --g 30 --f 1 --p 1

#probando generacion sintetica
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 27 34 --g 30 --f 1
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 0 --f 1
