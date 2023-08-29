import torch
import pandas as pd
import numpy as np

from SegmentationUNet import SegmentationUNet
from SegmentationDataset import SegmentationDataset
from torch.utils.data import DataLoader	

import os
import mrcfile

num_classes=3
folds = 3
depth=4
img_size=(64,64,64)
device = 'cuda' if torch.cuda.is_available() else 'cpu'
extra_width = 8
seed = 42
randg = torch.Generator()
randg.manual_seed(seed)


df = pd.read_csv('dataset/dataset_patches_final_new.csv', dtype=str)
test_list = ['4671', '5140', '6618', '8723', '22423', '8722', '4156', '4400', '20495', '13387', '20613', '11838', '10274', '24578', '25971', '3435', '8438', '25980', '3353', '6207', '4057', '11688', '22964', '3355', '30239', '22754', '21693', '22115', '25972', '23064', '8721']
# Get list of entry_id, subunit 
unique_df = df.drop_duplicates(subset=['id', 'subunit'])
grouped_df = df.groupby(['id','subunit'])
test_split = df[df['id'].isin(test_list)]
test_split.reset_index(inplace=True)

test_dataset = SegmentationDataset(test_split, num_classes, img_size, randg, device, augmentate=False, extra_width=extra_width, is_validation=True, is_nonoverlap_stride=True)

weights = 'results/3DUnet-0x8fb3-_fold-0-1_20230606-155914/best_model_198_0.8113.pt'
test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)

model = SegmentationUNet(num_classes, device, in_channels=2, depth=depth, batch_norm=True)
model.load_state_dict(torch.load(weights, map_location=torch.device(device)))
model.to(device)
model.eval()

for inputs_y,row_ix in zip(test_loader,test_split.index.tolist()):
    row = test_split.iloc[[row_ix]]
    inputs,y = inputs_y
    with torch.no_grad():
        pred = model(inputs)
    prob = torch.nn.functional.softmax(pred, dim=1)
    print("Saving {}_{}_{}".format(row['id'].item(),row['subunit'].item(),row['patch'].item()))
    inputs = inputs.detach().cpu().numpy()
    densities = inputs[:,0,:]
    points = inputs[:,1,:]
    map_id = row['id'].item()
    segment_id = row['subunit'].item()
    patch_ix = row['patch'].item()
    #np.save("out_new/map_{}_{}_{}.npy".format(map_id, segment_id, patch_ix),densities)
    #np.save("out_new/points_{}_{}_{}.npy".format(map_id, segment_id, patch_ix),points)
    #np.save("out_new/gt_{}_{}_{}.npy".format(map_id, segment_id, patch_ix),y.detach().cpu().numpy())
    np.save("out_new/logits_{}_{}_{}.npy".format(map_id, segment_id, patch_ix),pred[0][2]) 
    preds = torch.argmax(prob, dim=1)
    #np.save("out_new/pred_{}_{}_{}.npy".format(map_id, segment_id, patch_ix),preds.detach().cpu().numpy())
    '''
    with mrcfile.open('out_new/pred_{}_{}_{}.mrc'.format(map_id, segment_id, patch_ix), mode='w+') as mrc:
        mrc.set_data(np.float32(preds.detach().cpu().numpy()[0]))
    with mrcfile.open('out_new/map_{}_{}_{}.mrc'.format(map_id,segment_id, patch_ix), mode='w+') as mrc:
        mrc.set_data(np.float32(densities[0]))
    with mrcfile.open('out_new/points_{}_{}_{}.mrc'.format(map_id,segment_id, patch_ix), mode='w+') as mrc:
        mrc.set_data(np.float32(points[0]))
    with mrcfile.open('out_new/gt_{}_{}_{}.mrc'.format(map_id,segment_id, patch_ix), mode='w+') as mrc:
        mrc.set_data(np.float32(y.detach().cpu().numpy()[0]))
    '''


for name, group in test_split.groupby(['id','subunit']):
    num_patches = len(group)
    map_id = name[0]
    segment_id = name[1]
    print("Loading patches for {} {}".format(map_id,segment_id))
    np_sample = np.load('/work/mzumbado/data_patches/test/{}_{}.npy'.format(map_id,segment_id))
    sample = torch.from_numpy(np_sample)[0,:,:,:][None,:]    
    subunit_size = np_sample.shape
    patches = sample.unfold(1, 64, 64).unfold(2, 64, 64).unfold(3, 64, 64)
    unfold_shape = patches.size()
    patches = patches.contiguous().view(-1,patches.size(0), 64, 64, 64)
    patches_pred = [] 
    patches_points = []
    patches_gt = []
    patches_map = []
    patches_logit = []
    for n in range(num_patches):
        patch_pred = np.load('out_new/pred_{}_{}_{}.npy'.format(map_id,segment_id,n))
        patch_map = np.load('out_new/map_{}_{}_{}.npy'.format(map_id,segment_id,n))
        patch_points = np.load('out_new/points_{}_{}_{}.npy'.format(map_id,segment_id,n))
        patch_gt = np.load('out_new/gt_{}_{}_{}.npy'.format(map_id,segment_id,n))
        patch_logit = np.load('out_new/logits_{}_{}_{}.npy'.format(map_id,segment_id,n))
        patches_pred.append(torch.from_numpy(patch_pred))
        patches_points.append(torch.from_numpy(patch_points))
        patches_gt.append(torch.from_numpy(patch_gt))
        patches_map.append(torch.from_numpy(patch_map))
        patches_logit.append(torch.from_numpy(patch_logit))
    patches_pred = torch.stack(patches_pred)
    patches_points = torch.stack(patches_points)
    patches_gt = torch.stack(patches_gt)
    patches_map=torch.stack(patches_map)
    patches_logit = torch.stack(patches_logit)
    #patches_reconstruct = torch.vstack([patches_gt,patches_map,patches_points,patches_pred])
    # reshape output to match F.fold input
    patches_pred = patches_pred.contiguous().view(-1, 64, 64, 64)
    patches_points = patches_points.contiguous().view(-1, 64, 64, 64)
    patches_gt = patches_gt.contiguous().view(-1, 64, 64, 64)
    patches_map = patches_map.contiguous().view(-1, 64, 64, 64) 
    patches_logit = patches_logit.contiguous().view(-1, 64, 64, 64)

    # Reshape back
    output_c = unfold_shape[1] * unfold_shape[4] 
    output_h = unfold_shape[2] * unfold_shape[5] 
    output_w = unfold_shape[3] * unfold_shape[6] 

    assert output_c == subunit_size[1]
    assert output_h == subunit_size[2]
    assert output_w == subunit_size[3]

    patches_pred = patches_pred.view(unfold_shape)
    patches_pred = patches_pred.permute(0, 1, 4, 2, 5, 3, 6).contiguous()
    patches_pred = patches_pred.view(1, output_c, output_h, output_w)
    patches_points = patches_points.view(unfold_shape)
    patches_points = patches_points.permute(0, 1, 4, 2, 5, 3, 6).contiguous()
    patches_points = patches_points.view(1, output_c, output_h, output_w)
    patches_map = patches_map.view(unfold_shape)
    patches_map = patches_map.permute(0, 1, 4, 2, 5, 3, 6).contiguous()
    patches_map = patches_map.view(1, output_c, output_h, output_w)
    patches_gt = patches_gt.view(unfold_shape)
    patches_gt = patches_gt.permute(0, 1, 4, 2, 5, 3, 6).contiguous()
    patches_gt = patches_gt.view(1, output_c, output_h, output_w)
    patches_logit = patches_logit.view(unfold_shape)
    patches_logit = patches_logit.permute(0, 1, 4, 2, 5, 3, 6).contiguous()
    patches_logit = patches_logit.view(1, output_c, output_h, output_w)
    #with mrcfile.open('out_new/full/pred_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
    #    mrc.set_data(np.float32(patches_pred))
    with mrcfile.open('out_new/full/logits_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
        mrc.set_data(np.float32(patches_logit))
    '''
    with mrcfile.open('out_new/full/map_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
        mrc.set_data(np.float32(patches_map))
    with mrcfile.open('out_new/full/points_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
        mrc.set_data(np.float32(patches_points))
    with mrcfile.open('out_new/full/gt_{}_{}.mrc'.format(map_id,segment_id), mode='w+') as mrc:
        mrc.set_data(np.float32(patches_gt))
    '''
