import numpy as np 
import torch
from torch.utils.data import Dataset
from torch import from_numpy
import pandas as pd

        
        
         

class SegmentationDataset(Dataset):
    def __init__(self, df, num_classes, image_size, randg, device,augmentate=False, extra_width=0, is_validation=False, is_nonoverlap_stride=False):
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
        self.is_nonoverlap_stride = is_nonoverlap_stride
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
        #print("Getting output for {} {} {}".format(map_id,segment_id,patch_ix))
         
        if self.is_nonoverlap_stride:
            stride = self.image_size[0]
            #datapath ='/work/mzumbado/data_patches/test/{}_{}.npy'.format(map_id,segment_id)
            datapath ='/work/mzumbado/patches/test/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix)
        else:
            stride = self.image_size[0] //2
            #datapath ='/work/mzumbado/data_patches/{}_{}.npy'.format(map_id,segment_id)
            datapath = '/work/mzumbado/patches/{}_{}_{}.torch'.format(map_id,segment_id,patch_ix)

        #np_array = np.load(datapath)
        #all_tensor = torch.from_numpy(np_array)
        all_tensor = torch.load(datapath)
        mask_data = all_tensor[0,:]
        map_data =  all_tensor[1,:]
        points_data = all_tensor[2:,:]
        
        # Load points according pool sample size
        if self.is_validation:
            point_idx = -1
        else:
            point_idx = np.random.choice(a=points_data.size(0), size=1 ).item()
        
        point_data = points_data[point_idx]/3
        '''
        mask_patches = mask_data.unfold(0, self.image_size[0], stride).unfold(1,self.image_size[0],stride).unfold(2,self.image_size[0],stride) 
        mask_patches = mask_patches.contiguous().view(-1,self.image_size[0],self.image_size[1],self.image_size[2])
        map_patches = map_data.unfold(0, self.image_size[0], stride).unfold(1,self.image_size[0],stride).unfold(2,self.image_size[0],stride)    
        map_patches = map_patches.contiguous().view(-1,self.image_size[0],self.image_size[1],self.image_size[2])
        point_patches = point_data.unfold(0, self.image_size[0], stride).unfold(1,self.image_size[0],stride).unfold(2,self.image_size[0],stride)    
        point_patches = point_patches.contiguous().view(-1,self.image_size[0],self.image_size[1],self.image_size[2])

        #print("File shape {}, folded {}".format(np_array.shape, patches_tensor.shape))
        mask_patch = mask_patches[patch_ix]
        map_patch = map_patches[patch_ix]
        point_patch = point_patches[patch_ix]
        '''
        mask_patch = mask_data
        map_patch = map_data
        point_patch = point_data
        
        try:
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
            # Create two channel input data
            input_data = torch.stack([map_patch, point_patch])
            #input_data = torch.stack([map_patch])
            x = input_data.to(dtype=torch.float)
            #mask_foreground = mask_patch==2.0
            #y = mask_foreground.to(dtype=torch.long)
            y = mask_patch.to(dtype=torch.long)
            del map_patch
            del mask_patch
            del point_patch
            #del mask_patches
            #del map_patches
            #del point_patches
            del all_tensor
            #del np_array
            if self.augmentate:
                x,y = self.transform(x,y)
            return x,y
        except Exception as e:
            print(e)
            

