#!/bin/bash
#SBATCH --job-name protein
#SBATCH --time=48:00:00
#SBATCH --partition nukwa
#SBATCH --nodelist=nukwa-05.cnca
#SBATCH --output=%j.%N.out


source ~/.bashrc
module load cuda/10.1.105
module load cudnn-10/10.1
export PYTHONPATH=$PYTHONPATH:/work/mzumbado/biostructure
conda activate /work/mzumbado/EMAN2/envs/pytorch

python run.py 
