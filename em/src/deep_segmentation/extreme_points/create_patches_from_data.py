import os
import torch
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F
import pandas as pd
import mrcfile



def extract_patches(test_id, segment_id):
    output_df = pd.DataFrame(columns=['id','subunit','patch','data_path']) 
    data_row = df_shape[(df_shape['id']==test_id) & (df_shape['chain_label']==segment_id+'.0')]
    map_shape = data_row['map_length'].item()
    map_shape = tuple([int(token) for token in map_shape.split(',')])

    row = unique_df[(unique_df['id']==test_id) & (unique_df['subunit']==segment_id)][['min_x','min_y','min_z','max_x','max_y','max_z']]


    min_x = max(int(row['min_x'])-extra_width, 0 )
    min_y = max(int(row['min_y'])-extra_width, 0)
    min_z = max(int(row['min_z'])-extra_width, 0)
    max_x = min(int(row['max_x'])+extra_width, map_shape[0])
    max_y = min(int(row['max_y'])+extra_width, map_shape[1])
    max_z = min(int(row['max_z'])+extra_width, map_shape[2])

    subunit_shape = tuple([max_x-min_x, max_y-min_y, max_z-min_z])

    # Compute original pad, slice and create corrected padding
    pad0_left = (subunit_shape[0] // stride_train * stride_train + size) - subunit_shape[0]
    pad1_left = (subunit_shape[1] // stride_train * stride_train + size) - subunit_shape[1]
    pad2_left = (subunit_shape[2] // stride_train * stride_train + size) - subunit_shape[2]
    # Calculate symmetric padding
    pad0_right = pad0_left // 2 if pad0_left % 2 ==0 else pad0_left // 2 + 1
    pad1_right = pad1_left // 2 if pad1_left % 2 ==0 else pad1_left // 2 + 1
    pad2_right = pad2_left // 2 if pad2_left % 2 ==0 else pad2_left // 2 + 1
    pad0_left = pad0_left // 2
    pad1_left = pad1_left // 2
    pad2_left = pad2_left // 2

    #remove original padding by slicing 

    array = np.load('/work/mzumbado/data_patches/{}_{}.npy'.format(test_id,segment_id))
    unpadded = array[:,pad0_left:array.shape[0]-pad0_right,pad1_left:array.shape[1]-pad1_right,pad2_left:array.shape[2]-pad2_right]


    # New padding
    pad0_left = (unpadded.shape[1] // stride_val * stride_val + size) - unpadded.shape[1]
    pad1_left = (unpadded.shape[2] // stride_val * stride_val + size) - unpadded.shape[2]
    pad2_left = (unpadded.shape[3] // stride_val * stride_val + size) - unpadded.shape[3]
    # Calculate symmetric padding
    pad0_right = pad0_left // 2 if pad0_left % 2 ==0 else pad0_left // 2 + 1
    pad1_right = pad1_left // 2 if pad1_left % 2 ==0 else pad1_left // 2 + 1
    pad2_right = pad2_left // 2 if pad2_left % 2 ==0 else pad2_left // 2 + 1
    pad0_left = pad0_left // 2
    pad1_left = pad1_left // 2
    pad2_left = pad2_left // 2

    
    padded = np.pad(unpadded, ((0,0),(pad0_left, pad0_right), (pad1_left, pad1_right), (pad2_left, pad2_right)))
   
    print("Saving data_patches/test/{}_{}.npy".format(test_id,segment_id))
    print("Shape: New {}, Old{}".format(padded.shape,array.shape))
    np.save('/work/mzumbado/data_patches/test/{}_{}.npy'.format(test_id,segment_id), padded)

    patches = torch.from_numpy(padded).unfold(3, size, stride_val).unfold(2, size, stride_val).unfold(1, size, stride_val) 

    patches = patches.contiguous().view(-1,12,size,size,size)
    number_patches = patches.size(0)
    print(patches.shape)

    for i in range(number_patches):
        output_df = output_df.append({'id':map_id, 'subunit':segment_id,'patch':i,'data_path':'data_patches/test/{}_{}.npy'.format(test_id,segment_id)}, ignore_index=True)
    # Return dataframe of created patches
    return output_df

    


df = pd.read_csv('dataset/dataset_extreme_points.csv',dtype=str)
df_shape = pd.read_csv('dataset/dataset_exp_merged.csv', dtype=str)
unique_df = df.drop_duplicates(subset=['id', 'subunit'])

output_df = pd.DataFrame(columns=['id','subunit','patch','data_path'])

stride_train = 32
stride_val = 64
size = 64
extra_width=8
output_path= 'dataset/test_data/'

df_patches = pd.read_csv('dataset/dataset_patches_final_new.csv', dtype=str)
test_list = ['4671', '5140', '6618', '8723', '22423', '8722', '4156', '4400', '20495', '13387', '20613', '11838', '10274', '24578', '25971', '3435', '8438', '25980', '3353', '6207', '4057', '11688', '22964', '3355', '30239', '22754', '21693', '22115', '25972', '23064', '8721']
#test_list = ['11688']
# Get list of entry_id, subunit 
grouped_df = df.groupby(['id','subunit'])
test_split = df_patches[df_patches['id'].isin(test_list)]
test_split.reset_index(inplace=True)

'''
for name, group in test_split.groupby(['id','subunit']):
    print(name)
    map_id = name[0]
    segment_id = name[1]
    output_df = output_df.append(extract_patches(map_id, segment_id), ignore_index=True)

output_df.to_csv('dataset/dataset_patches_test.csv', index=False)
'''
train_val_split = df_patches[~df_patches['id'].isin(test_list)]
train_val_split.reset_index(inplace=True)
'''
final_patches_df = train_val_split.append(output_df, ignore_index=True)
final_patches_df.to_csv('dataset/dataset_patches_final_new.csv', index=False)
'''

# Get train stats

for i,row in train_val_split.iterrows():
    map_id = row[1]
    segment_id = row[2]
    patch_ix = int(row[3])
    if os.path.exists('/work/mzumbado/patches/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix)):
        continue
    print("Loading training EMDB-{}, subunit {}, patch {}".format(map_id,segment_id,patch_ix))
    array = torch.from_numpy(np.load('/work/mzumbado/data_patches/{}_{}.npy'.format(map_id,segment_id)))
    mask_data = array[0,:]
    map_data = array[1,:]
    point_data = array[2:11,:]

    mask_patches = mask_data.unfold(0, size, stride_train).unfold(1,size,stride_train).unfold(2,size,stride_train)
    mask_patches = mask_patches.contiguous().view(-1,size,size,size)
    map_patches = map_data.unfold(0, size, stride_train).unfold(1,size,stride_train).unfold(2,size,stride_train)
    map_patches = map_patches.contiguous().view(-1,size,size,size)
    point_patches = point_data.unfold(1, size, stride_train).unfold(2,size,stride_train).unfold(3,size,stride_train)
    point_patches = point_patches.contiguous().view(9,-1,size,size,size)

    mask_patch = mask_patches[patch_ix]
    map_patch = map_patches[patch_ix]
    point_patch = point_patches[:,patch_ix]

    stack_list = [mask_patch, map_patch]
    point_patches =  [ p for p in point_patch] 
    stack_list = stack_list + point_patches
    data = torch.stack(stack_list)

    torch.save(data, '/work/mzumbado/patches/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix))

for i,row in test_split.iterrows():
    map_id = row[1]
    segment_id = row[2]
    patch_ix = int(row[3])
    if os.path.exists('/work/mzumbado/patches/test/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix)):
        continue
    print("Loading test EMDB-{}, subunit {}, patch {}".format(map_id,segment_id,patch_ix))
    array = torch.from_numpy(np.load('/work/mzumbado/data_patches/test/{}_{}.npy'.format(map_id,segment_id)))
    mask_data = array[0,:]
    map_data = array[1,:]
    point_data = array[2:11,:]

    mask_patches = mask_data.unfold(0, size, stride_val).unfold(1,size,stride_val).unfold(2,size,stride_val)
    mask_patches = mask_patches.contiguous().view(-1,size,size,size)
    map_patches = map_data.unfold(0, size, stride_val).unfold(1,size,stride_val).unfold(2,size,stride_val)
    map_patches = map_patches.contiguous().view(-1,size,size,size)
    point_patches = point_data.unfold(1, size, stride_val).unfold(2,size,stride_val).unfold(3,size,stride_val)
    point_patches = point_patches.contiguous().view(9,-1,size,size,size)

    mask_patch = mask_patches[patch_ix]
    map_patch = map_patches[patch_ix]
    point_patch = point_patches[:,patch_ix]

    stack_list = [mask_patch, map_patch]
    point_patches =  [ p for p in point_patch]
    stack_list = stack_list + point_patches
    data = torch.stack(stack_list)

    torch.save(data, '/work/mzumbado/patches/test/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix))

