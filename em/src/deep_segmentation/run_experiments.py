
import numpy as np
import os


nodes_v = [ 'nukwa-04.cnca', 'nukwa-05.cnca', 'nukwa-06.cnca']
nodes_k = [ 'nukwa-00.cnca', 'nukwa-01.cnca', 'nukwa-02.cnca' , 'nukwa-03.cnca']

sizes = [ 72, 80, 88, 96, 104, 112, 120, 128 ]
losses = ['CrossEntropy', 'Weigthed', 'Tversky', 'Dice']
batch_sizes = [ 1, 2 ]
optims = ['SGD', 'Adam']
lrs = np.logspace(-4, -2, num=5)
decays = np.logspace(-6, -4, num=5) 
widths = [0,4,8,12]
alphas = [0.5, 0.6, 0.7, 0,8, 0.9]
for j in range(3):
    base_content = """#!/bin/bash
#SBATCH --job-name deep_protein 
#SBATCH --time=7-00:00:00
#SBATCH --partition nukwa-long
#SBATCH --output=%j.%N.out
##SBATCH --nodes=1
##SBATCH --tasks-per-node=1
##SBATCH --ntasks=1

source ~/.bashrc
module load cuda/10.1.105
module load cudnn-10/10.1
export MKL_THREADING_LAYER=GNU
export PYTHONPATH=$PYTHONPATH:/work/mzumbado/biostructure/em/src
{}


"""
    with open('train.slurm', 'w') as jobfile:
        size = np.random.choice(sizes, size=1)[0]
        loss = np.random.choice(losses, size=1)[0]
        lr = np.random.choice(lrs, size=1)[0]
        optim = np.random.choice(optims, size=1)[0]
        decay = np.random.choice(decays, size=1)[0]
        width = np.random.choice(widths, size=1)[0]
        batch_size = np.random.choice(batch_sizes, size=1)[0]
        alpha = np.random.choice(alphas, size=1)[0]
        beta = 1 - alpha
        if (size > 72 | (size == 72 & batch_size == 2)):
            node = np.random.choice(nodes_v, size=1)[0]
            base_content = base_content.format('conda activate /work/mzumbado/EMAN2/envs/pytorch\n')
        else:
            node = np.random.choice(nodes_k, size=1)[0]
            base_content = base_content.format('conda activate /work/mzumbado/EMAN2/envs/pytorch_kepler\n')
        if loss in ['Weigthed', 'Tversky']:
            cmd = "python run_.py run --with_amp=False --input_size={} --batch_size={} --alpha={} --beta={} --num_epochs=150 --optim={} --learning_rate={} --weight_decay={} --criterion={} --extra_width={} --validate_every=1".format(size,batch_size,alpha,beta,optim,lr,decay,loss,width)
        else:
            cmd = "python run_.py run --with_amp=False --input_size={} --batch_size={} --alpha={} --beta={} --num_epochs=100 --optim={} --learning_rate={} --weight_decay={} --criterion={} --extra_width={} --validate_every=1".format(size,batch_size,alpha,beta,optim,lr,decay,loss,width)
        print(base_content)
        print(cmd)
        n =jobfile.write(base_content)
        n =jobfile.write(cmd)
    cmd = "sbatch --nodelist={} train.slurm".format(node)
    print(cmd)
    os.system(cmd)

