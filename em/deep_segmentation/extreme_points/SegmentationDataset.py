import numpy as np 
from torch.utils.data import Dataset
from torch import from_numpy

from scipy.ndimage import zoom
import pandas as pd
from em.molecule import Molecule

class SegmentationDataset(Dataset):
    def __init__(self, df, num_classes, image_size, device):
        """
        Dataset class for EM data 
        :param num_classes: number of classes to classify
        :param device: CPU or GPU
        """
        self.maps = df['map_path'].tolist() 
        self.contours = df['contourLevel'].tolist()
        self.points = df['tagged_points_path'].tolist()
        self.masks = df['tagged_path'].tolist()
        self.num_classes = num_classes
        self.image_size = image_size
        self.device = device

    def __len__(self):
        return len(self.maps)

    def __getitem__(self, idx):
        map_data = Molecule(self.maps[idx], self.contours[idx], [1.2]).getDataAtContour(1.2)
        mask_data = np.load(self.masks[idx])
        point_data =  np.load(self.points[idx])
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
        norm_data = (map_data - data_min)/ (data_max-data_min)
        # Create two channel input data
        input_data = np.vstack([norm_data[np.newaxis], point_data[np.newaxis]])

        x = from_numpy(input_data).float().to(self.device)
        y = from_numpy(mask_data).long().to(self.device)
        return x, y
