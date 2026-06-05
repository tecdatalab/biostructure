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
#SBATCH --tasks-per-node=30
#SBATCH --partition=dribe
#SBATCH --ntasks=30
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=dribe-01.cnca

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
##conda init bash 

# Load required modules
module load gcc/9.3.0
module load mpich/3.3.2-gcc-9.3.0

# Setup conda
eval "$(/home/ggutierrez/miniconda/bin/conda shell.bash hook)"
source ~/.bashrc

# Create and activate environment
conda create -n visual python=3.6 -y
conda activate visual

# Install required packages
conda install -y -c conda-forge numpy pandas mpi4py mrcfile scipy scikit-image pytorch

# Set Python path
export PYTHONPATH=$PYTHONPATH:/work/ggutierrez/biostructure

# Run feature extraction with specific result directory
time srun --mpi=pmi2 --cpu-bind=cores python -m mpi4py.futures feature_extraction.py --p 3 --n 3 --result_dir /data/ggutierrez/extreme_points
