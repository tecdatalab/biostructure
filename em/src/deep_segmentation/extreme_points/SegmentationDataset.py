import numpy as np 
import torch
from torch.utils.data import Dataset
from torch import from_numpy
import pandas as pd
from scipy.ndimage import zoom
import mrcfile

        
        
         

class SegmentationDataset(Dataset):
    def __init__(self, df, num_classes, image_size, randg, device,augmentate=False, extra_width=0, is_validation=False):
        """
        Dataset class for EM data 
        :param num_classes: number of classes to classify
        """
        self.unique_dataframe = df.drop_duplicates(subset=['id','subunit']).reset_index(drop=True)
        self.points_df = df
        self.maps = self.unique_dataframe['map_path'].tolist() 
        self.contours = self.unique_dataframe['contourLevel'].tolist()
        self.masks = self.unique_dataframe['tagged_path'].tolist()
        self.extra_width = extra_width
        self.bbox_coords=  (self.unique_dataframe['min_x'].tolist(),self.unique_dataframe['min_y'].tolist(),self.unique_dataframe['min_z'].tolist(),self.unique_dataframe['max_x'].tolist(),self.unique_dataframe['max_y'].tolist(),self.unique_dataframe['max_z'].tolist())
        self.num_classes = num_classes
        self.image_size = image_size
        self.randg = randg
        self.augmentate = augmentate
        self.is_validation = is_validation
        self.device = device

    def __len__(self):
        return len(self.unique_dataframe)

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
        map_id = self.unique_dataframe.loc[[idx]]['id'].item()
        segment_id = self.unique_dataframe.loc[[idx]]['subunit'].item()
        map_object = mrcfile.open(self.maps[idx])
        map_data = np.copy(map_object.data)
        map_data[map_data<self.contours[idx]]=0
        mask_data = np.load(self.masks[idx]) 
        min_x = max(self.bbox_coords[0][idx]-self.extra_width, 0 )
        min_y = max(self.bbox_coords[1][idx]-self.extra_width, 0)
        min_z = max(self.bbox_coords[2][idx]-self.extra_width, 0)
        max_x = min(self.bbox_coords[3][idx]+self.extra_width, map_data.shape[0])
        max_y = min(self.bbox_coords[4][idx]+self.extra_width, map_data.shape[1])
        max_z = min(self.bbox_coords[5][idx]+self.extra_width, map_data.shape[2])
        try:
            # Remove noise
            #map_data[mask_data==0] = 0
            mask_data = mask_data[min_x:max_x,min_y:max_y,min_z:max_z]
            map_data = map_data[min_x:max_x,min_y:max_y,min_z:max_z]
            # Load points according pool sample size
            if self.is_validation:
                points_df = self.points_df[(self.points_df['id']==map_id) & (self.points_df['subunit']==segment_id)]
                point_data_filename = points_df['tagged_points_path'].iloc[-1]
            else:
                points_df = self.points_df[(self.points_df['id']==map_id) & (self.points_df['subunit']==segment_id)][:-1]
                point_data_filename = points_df.sample(1, random_state=self.randg.seed() % 2**(32)-1)['tagged_points_path'].item()
            #print("fetching map {} subunit {} points {}".format(map_id,segment_id,point_data_filename))
            point_data =  np.load(point_data_filename)
            point_data = point_data[min_x:max_x,min_y:max_y,min_z:max_z]
            # Resize imput data
            # zoom_factor = [ resized_shape/axis_shape for axis_shape,resized_shape in zip(map_data.shape,self.image_size) ]
            #map_data = zoom(map_data, zoom_factor, order=1)
            # Nearest neighbor interpolation
            #mask_data = zoom(mask_data, zoom_factor, order=0)
            # Nearest neighbor interpolation
            #point_data = zoom(point_data, zoom_factor, order=0)
            # Input data normalization
            data_max = np.max(map_data)
            data_min = np.min(map_data)
            norm_data = (map_data - data_min)/ (data_max-data_min + 1e-6)
            # Create two channel input data
            input_data = np.vstack([norm_data[np.newaxis], point_data[np.newaxis]])
            '''
            # Zero pad to match shape
            pad = [[0,0],[0,0], [0,0], [0,0]] 
            pad_zeros = [ s - self.image_size[0] for s in map_data.shape]
            #Check padding
            for i,c in enumerate(pad_zeros):
                if c <= 0:
                    pad_len = self.image_size[0] - norm_data.shape[i]
                    pad[i+1][1] = pad_len
                    pad_zeros[i] = 0
            if np.sum(pad)>0:
                input_data = np.pad(input_data, pad)
                mask_data = np.pad(mask_data, pad[1:])
            '''
            pad_zeros = [ (self.image_size[0]-s) % self.image_size[0] for s in map_data.shape]
            even_pad = [True if p%2==0 else False for p in pad_zeros]
            pad = [ [p//2,p//2] if is_even else [p//2,p//2+1] for p,is_even in zip(pad_zeros,even_pad) ] 
            pad= [[0,0]]+pad
            input_data = np.pad(input_data,pad)
            mask_data = np.pad(mask_data,pad[1:]) 
            x = from_numpy(input_data).float()
            y = from_numpy(mask_data).long()
            if self.augmentate:
                x,y = self.transform(x,y)
            x_patches = x.unfold(3, self.image_size[2], 48).unfold(2, self.image_size[1], 48).unfold(1, self.image_size[0], 48)
            x_patches = x_patches.contiguous().view(-1, 2, 32, 32, 32)
            y_patches = y.unfold(2, 96, 48).unfold(1, 96, 48).unfold(0, 96, 48)
            y_patches = y_patches.contiguous().view(-1, 32, 32, 32)
            indices = torch.randperm(x_patches.shape[0])
            x_patches=x_patches[indices]
            y_patches=y_patches[indices]
            return x_patches,y_patches
        except Exception as e:
            print(e)
            print("dim0 [{},{}], dim1 [{},{}], dim2 [{},{}]".format(min_x,max_x,min_y,max_y,min_z,max_z))
            print(map_data.shape)
            print(self.image_size)
        map_object.close()
        del map_data
        del mask_data

