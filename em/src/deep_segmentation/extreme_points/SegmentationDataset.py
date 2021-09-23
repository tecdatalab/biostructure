import numpy as np 
import torch
from torch.utils.data import Dataset
from torch import from_numpy

from scipy.ndimage import zoom, distance_transform_edt, gaussian_filter
import pandas as pd
from em.molecule import Molecule

        
        
         

class SegmentationDataset(Dataset):
    def __init__(self, df, num_classes, image_size, randg, extra_width=0):
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
        map_data = Molecule(self.maps[idx], self.contours[idx]).getDataAtContour(1)
        mask_data = np.load(self.masks[idx]) 
        min_x = max(self.bbox_coords[0][idx]-self.extra_width, 0 )
        min_y = max(self.bbox_coords[1][idx]-self.extra_width, 0)
        min_z = max(self.bbox_coords[2][idx]-self.extra_width, 0)
        max_x = min(self.bbox_coords[3][idx]+self.extra_width, map_data.shape[0]-1)
        max_y = min(self.bbox_coords[4][idx]+self.extra_width, map_data.shape[1]-1)
        max_z = min(self.bbox_coords[5][idx]+self.extra_width, map_data.shape[2]-1)
        try:
            # Remove noise
            map_data[mask_data==0] = 0
            mask_data = mask_data[min_x:max_x,min_y:max_y,min_z:max_z]
            map_data = map_data[min_x:max_x,min_y:max_y,min_z:max_z]
            # Load points according pool sample size
            points_df = self.points_df[(self.points_df['id']==map_id) & (self.points_df['subunit']==segment_id)]
            point_data_filename = points_df.sample(1)['tagged_points_path'].item()
            #print("fetching map {} subunit {} points {}".format(map_id,segment_id,point_data_filename))
            point_data =  np.load(point_data_filename)
            point_data = point_data[min_x:max_x,min_y:max_y,min_z:max_z]
            #point_data =  compute_points(mask_data, 3)
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
            # Create two channel input data
            input_data = np.vstack([norm_data[np.newaxis], point_data[np.newaxis]])
            x = from_numpy(input_data).float()
            y = from_numpy(mask_data).long()
            x,y = self.transform(x,y)
            return x,y
        except Exception as e:
            print(e)
            print("dim0 [{},{}], dim1 [{},{}], dim2 [{},{}]".format(min_x,max_x,min_y,max_y,min_z,max_z))
            print(map_data.shape)
            print(self.image_size)

def compute_points(mask_map, number_points=3, gaussian_std=3):
    #print("unique",np.unique(tagged_map))
        #print("pathh {}".format(region_path))
    distance = distance_transform_edt(mask_map)
    distance[distance != 1] = 0
    index_x, index_y, index_z = np.where(distance == 1)
    chosen_indexes = np.random.choice(len(index_x), number_points, replace=False)
    index_x = index_x[chosen_indexes]
    index_y = index_y[chosen_indexes]
    index_z = index_z[chosen_indexes]
    point_array = np.zeros_like(mask_map)
    point_array[index_x,index_y,index_z] = 1.0
    point_array = gaussian_filter(point_array, gaussian_std)
    return point_array 
