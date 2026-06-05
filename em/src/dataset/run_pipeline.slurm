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
#SBATCH --job-name=biostructure
#SBATCH --output=%j.%N.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=36
#SBATCH --tasks-per-node=36
#SBATCH --mem=1000GB
#SBATCH --partition=dribe-long
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gergutierrez@estudiantec.cr

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
##SBATCH --partition=nukwa-debug
##SBATCH --ntasks=48    
##SBATCH --time=04:00:00
##SBATCH --nodelist=nukwa-07.cnca

# Load required modules
module load gcc/9.3.0

# Check initial system state
echo "Job started at: $(date)"
echo "Checking initial disk space..."
df -h /tmp/$USER/conda /data/$USER || true
echo "Checking initial memory state..."
free -h

# Set temporary directories for conda
export CONDA_PKGS_DIRS=/tmp/$USER/conda/pkgs
export CONDA_ENVS_PATH=/tmp/$USER/conda/envs
mkdir -p $CONDA_PKGS_DIRS $CONDA_ENVS_PATH

# Clean existing environment and caches
rm -rf /tmp/$USER/conda/envs/biostructure
conda clean --all -y
rm -rf ~/.cache/pip

# Setup conda environment
eval "$(/home/$USER/miniconda/bin/conda shell.bash hook)"
source ~/.bashrc

# Create and activate environment with Python 3.10
conda create --prefix /tmp/$USER/conda/envs/biostructure python=3.10 numpy pandas -y
source activate /tmp/$USER/conda/envs/biostructure

# Install packages with progress checks
echo "Installing mrcfile..."
conda install -y -c conda-forge mrcfile && echo "mrcfile installed"

echo "Installing scipy and scikit-image..."
conda install -y -c conda-forge scipy scikit-image && echo "scipy and scikit-image installed"

echo "Installing pytorch and biopython..."
conda install -y -c conda-forge pytorch biopython && echo "pytorch and biopython installed"

# Verify installations
python -c "import numpy; import pandas; import mrcfile; import scipy; import torch; import Bio; print('All packages verified')"

# Create and prepare output directories
mkdir -p /home/$USER/extreme_points
mkdir -p /data/$USER/extreme_points
chmod 777 /home/$USER/extreme_points
chmod 777 /data/$USER/extreme_points

# Run feature extraction with monitoring
echo "Starting feature extraction at $(date)"
python feature_extraction.py --result_dir /home/$USER/extreme_points --final_dir /data/$USER/extreme_points
EXIT_CODE=$?

# Check final system state
echo ""
echo "Checking final disk space..."
df -h /tmp/$USER/conda /data/$USER || true
echo "Checking final memory state..."
free -h
echo "Feature extraction exit code: $EXIT_CODE"

# Cleanup
conda deactivate
rm -rf /tmp/$USER/conda/pkgs/*

exit $EXIT_CODE