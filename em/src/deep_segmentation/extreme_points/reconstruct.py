import torch
import pandas as pd
import numpy as np


import os
import mrcfile

size=64
stride_train = 32
stride_test = 64

df = pd.read_csv('dataset/dataset_patches_final.csv', dtype=str)
test_list = ['11688']
#test_list = ['3329', '3354', '4039', '5248', '20916', '5247', '31118', '31823', '13371', '12312', '13378', '11688', '10272', '22965', '23930', '2752', '5002', '8231', '25972', '3948', '31028', '5691', '20932', '25980', '24712', '20930', '20665', '7305', '8230', '22115', '5246']
# Get list of entry_id, subunit 
test_split = df[df['id'].isin(test_list)]
test_split.reset_index(inplace=True)

patches_gt = []
patches_map = []
patches_points = []

for i,row in test_split.iterrows():
    map_id = row[1]
    segment_id = row[2]
    patch_ix = int(row[3])
    print("Loading patch {} for EMDB {} subunit {}".format(patch_ix,map_id,segment_id))
    array = torch.load('/work/mzumbado/patches/test/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix))
    patch_map = array[1,:]
    patch_points = array[2,:]
    patch_gt = array[0,:]
    patches_points.append(patch_points)
    patches_gt.append(patch_gt)
    patches_map.append(patch_map)
    
np_sample = np.load('/work/mzumbado/data_patches/test/{}_{}.npy'.format(map_id,segment_id))
sample = torch.from_numpy(np_sample)[0,:,:,:]
print("Sample shape {}".format(sample.shape))
unfold_sample = sample.unfold(0,size,stride_test).unfold(1,size,stride_test).unfold(2,size,stride_test)
unfold_shape = unfold_sample.size()
print("unfolded {}".format(unfold_shape))

patches_points = torch.stack(patches_points)
patches_gt = torch.stack(patches_gt)
patches_map=torch.stack(patches_map)
# reshape output to match F.fold input
patches_points = patches_points.contiguous().view(-1, 64, 64, 64)
patches_gt = patches_gt.contiguous().view(-1, 64, 64, 64)
patches_map = patches_map.contiguous().view(-1, 64, 64, 64) 

# Reshape back
output_c = unfold_shape[0] * unfold_shape[3] 
output_h = unfold_shape[1] * unfold_shape[4] 
output_w = unfold_shape[2] * unfold_shape[5] 


patches_points = patches_points.view(unfold_shape)
patches_points = patches_points.permute(0, 3, 1, 4, 2, 5).contiguous()
patches_points = patches_points.view(1, output_c, output_h, output_w)
patches_map = patches_map.view(unfold_shape)
patches_map = patches_map.permute(0, 3, 1, 4, 2, 5).contiguous()
patches_map = patches_map.view(1, output_c, output_h, output_w)
patches_gt = patches_gt.view(unfold_shape)
patches_gt = patches_gt.permute(0, 3, 1, 4, 2, 5).contiguous()
patches_gt = patches_gt.view(1, output_c, output_h, output_w)

with mrcfile.open('out_new/full/map_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
    mrc.set_data(np.float32(patches_map))
with mrcfile.open('out_new/full/points_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
    mrc.set_data(np.float32(patches_points))
with mrcfile.open('out_new/full/gt_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
    mrc.set_data(np.float32(patches_gt))
    
