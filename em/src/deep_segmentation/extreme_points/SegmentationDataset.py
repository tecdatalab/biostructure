import numpy as np 
import torch
from torch.utils.data import Dataset
from torch import from_numpy
import pandas as pd

        
        
         

class SegmentationDataset(Dataset):
    def __init__(self, df, num_classes, image_size, randg, device,augmentate=False, extra_width=0, is_validation=False):
        """
        Dataset class for EM data 
        :param num_classes: number of classes to classify
        """
        self.patches_df = df
        self.num_classes = num_classes
        self.image_size = image_size
        self.randg = randg
        self.augmentate = augmentate
        self.is_validation = is_validation
        self.device = device

 
    def __len__(self):
        return len(self.patches_df)


    def transform(self, x_in, y_in):
        x_out = x_in
        y_out = y_in
        prob = torch.rand(1, generator=self.randg)
        if prob >= 0.5:
            prob = torch.rand(1, generator=self.randg)
            if prob >= 0.5:
                x_out = torch.flip(x_out, [-3])
                y_out = torch.flip(y_out, [-3])
            prob = torch.rand(1, generator=self.randg)
            if prob >=0.5:
                x_out = torch.flip(x_out, [-2])
                y_out = torch.flip(y_out, [-2])
            prob = torch.rand(1, generator=self.randg)
            if prob >= 0.5:
                x_out = torch.flip(x_out, [-1])
                y_out = torch.flip(y_out, [-1])
        prob = torch.rand(1, generator=self.randg)
        if prob >= 0.5:
            k_l = [ 1, 2, 3 ]
            axis = torch.tensor([-1, -2, -3])
            angle = k_l[torch.randperm(len(k_l), generator=self.randg)[0]]
            axis = tuple(axis[torch.randperm(len(axis), generator=self.randg)[0:2]].tolist())
            x_out = torch.rot90(x_out,angle, axis)
            y_out = torch.rot90(y_out, angle, axis)
        return x_out, y_out

    def __getitem__(self, idx):
        #idx = 153
        row = self.patches_df.iloc[[idx]]
        map_id = row['id'].item()
        segment_id = row['subunit'].item()
        patch_ix = int(row['patch'].item()) 
        np_sample = np.load('/work/mzumbado/data_patches/{}_{}.npy'.format(map_id,segment_id))
        sample = torch.from_numpy(np_sample)
        patches_tensor = sample.unfold(3, self.image_size[2], self.image_size[2]//2).unfold(2, self.image_size[1], self.image_size[1]//2).unfold(1, self.image_size[0], self.image_size[0]//2)
        patches_tensor = patches_tensor.contiguous().view(-1,12,self.image_size[0],self.image_size[1],self.image_size[2])
        patch_data = patches_tensor[patch_ix]
        mask_data = patch_data[0]
        map_data = patch_data[1]
        points_data = patch_data[2:12]
        try:
            # Load points according pool sample size
            if self.is_validation:
                point_idx = -1
            else:
                point_idx = np.random.choice(a=points_data.size(0), size=1 ).item()
            ''' 
            # Resize imput data
            zoom_factor = [ resized_shape/axis_shape for axis_shape,resized_shape in zip(map_data.shape,self.image_size) ]
            map_data = zoom(map_data, zoom_factor, order=1)
            # Nearest neighbor interpolation
            mask_data = zoom(mask_data, zoom_factor, order=0)
            # Nearest neighbor interpolation
            point_data = zoom(point_data, zoom_factor, order=0)
            
            # Input data normalization
            data_max = np.max(map_data)
            data_min = np.min(map_data)
            norm_data = (map_data - data_min)/ (data_max-data_min + 1e-6)
            '''
            point_data = points_data[point_idx]/3
            # Create two channel input data
            input_data = torch.stack([map_data, point_data])
            x = input_data.float()
            y = mask_data.long()
            if self.augmentate:
                x,y = self.transform(x,y)
            return x,y
        except Exception as e:
            print(e)
            del patches_tensor

