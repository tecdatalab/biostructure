#!/bin/bash
#SBATCH --job-name 3DUnet
#SBATCH --time=7-00:00:00
#SBATCH --partition nukwa-long
#SBATCH --nodelist=nukwa-04.cnca
#SBATCH --output=%j.%N.out

source ~/.bashrc

conda activate pytorch_gpu

python run_.py run --with_amp=False --input_size=96  --num_epochs=200 --optim=SGD --learning_rate=0.001 --gamma=0.6 --depth=4 --batch_size=16 --batch_norm=False --criterion=Unified --extra_width=8 --validate_every=1
