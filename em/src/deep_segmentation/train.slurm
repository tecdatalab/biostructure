#!/bin/bash
#SBATCH --job-name 3DUnet
#SBATCH --time=7-00:00:00
#SBATCH --partition nukwa-long
#SBATCH --nodelist=nukwa-06.cnca
#SBATCH --output=%j.%N.out

source ~/.bashrc

conda activate pytorch_gpu

python run_.py run --with_amp=False --input_size=64 --gamma=0.5 --num_epochs=200 --optim=Adam --learning_rate=0.0008 --depth=4 --batch_size=64 --criterion=Dice  --extra_width=8 --validate_every=5 --weight_decay=0.00001 --seed=4325 --batch_norm=True --resume_from=results/3DUnet-0x89a1-_fold-2-1_20230721-130352/training_checkpoint_67600.pt
