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
#SBATCH --job-name=protein
#SBATCH --output=%j.%N.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=dribe
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=dribe-00.cnca

#ANDALAN
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=1
##SBATCH --tasks-per-node=48
##SBATCH --partition=andalan
##SBATCH --ntasks=48
##SBATCH --time=24:00:00

##NUKWA
##SBATCH --job-name=proteins
##SBATCH --output=%j.%N.out
##SBATCH --nodes=1
##SBATCH --tasks-per-node=48
##SBATCH --partition=nukwa
##SBATCH --ntasks=48    
##SBATCH --time=1-00:00:00
##SBATCH --nodelist=nukwa-04.cnca
#conda init bash 
source ~/.bashrc
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 1  --v 1 --s 22 29 --g 30 --f 1 --p 1

#probando generacion sintetica
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 1 --s 27 34 --g 30 --f 1
#time srun --mpi=pmi2 --cpu-bind=cores /home/mzumbado/.conda/envs/visual/bin/python3.6 -m  mpi4py.futures process_data.py --d 0 --v 0 --f 1
conda activate visual
export PYTHONPATH=$PYTHONPATH:/work/mzumbado/biostructure

time python main.py --input_path /work/mzumbado/biostructure/em/evolutionary_segmentation/*8528*.mrc --output_path /work/mzumbado/biostructure/em/dataset/results_3 --level 0.001
#time python segger_comparison.py
