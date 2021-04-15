import numpy as np
from torch.optim import Adam
from torch.utils.data import DataLoader

from SegmentationDataset import SegmentationDataset
from SegmentationUNet import SegmentationUNet
from CrossEntropyLoss import CrossEntropyLoss

from em.molecule import Molecule
from scipy.ndimage import zoom

import pandas as pd

def compute_class_weights(df, labels):
    """
    Helper method to compute class weights for training data
    :param df: Dataframe of samples
    :param label: List of numerical class labels to compute weights
    :return: tensor of weights for corresponding label
    """
    masks = df['tagged_path'].tolist()
    weights = np.zeros(len(labels))
    total = 0
    for mask in masks:
        mask_array = np.load(mask)
        w = [ np.sum(mask_array==l) for l in labels ]
        total += np.sum(w)
        weights += np.array(w)
    return total/weights
        
    

class SegmentationAgent:
    def __init__(self, val_percentage, test_num, num_classes,
                 batch_size, img_size, csv_path, shuffle_data,
                 learning_rate, depth, device):
        """
        A helper class to facilitate the training of the model
        """
        self.device = device
        self.num_classes = num_classes
        self.batch_size = batch_size
        self.img_size = img_size
        self.dataframe = pd.read_csv(csv_path)
        self.depth = depth
        train_split, val_split, test_split = self.make_splits(val_percentage, test_num, shuffle_data)
        #train_split, val_split = self.make_splits(val_percentage, test_num, shuffle_data)
        self.train_loader = self.get_dataloader(train_split)
        self.validation_loader = self.get_dataloader(val_split)
        self.test_loader = self.get_dataloader(test_split,'test')
        self.model = SegmentationUNet(self.num_classes, self.device, depth=self.depth)
        self.class_weights = compute_class_weights(train_split, [0,1])
        print("Class weights {}".format(self.class_weights))
        self.criterion = CrossEntropyLoss(self.num_classes, self.class_weights, self.device)
        self.optimizer = Adam(self.model.parameters(), lr=learning_rate)
        self.model.to(self.device)

    def make_splits(self, val_percentage=0.2, test_num=10, shuffle=True):
        """
        Split the data into train, validation and test datasets
        :param val_percentage: A decimal number which tells the percentage of
                data to use for validation
        :param test_num: The number of images to use for testing
        :param shuffle: Shuffle the data before making splits
        :return: tuples of splits
        """
        if shuffle:
            self.dataframe = self.dataframe.sample(frac=1)

        val_num = int(val_percentage * len(self.dataframe))
        print("val num {}".format(val_num))
        train_samples = self.dataframe.iloc[:val_num]
        validation_samples = self.dataframe.iloc[val_num:-test_num]
        test_samples = self.dataframe.iloc[-test_num:]
        print("Number of samples: {} Training {} Validation {} Testing {}".format((len(self.dataframe), len(train_samples), len(validation_samples), len(test_samples)))
        return train_samples, validation_samples, test_samples

    def get_dataloader(self, split, split_name=None):
        """
        Create a DataLoader for the given split
        :param split: train split, validation split or test split of the data
        :return: DataLoader
        """
        dataset = SegmentationDataset(split, self.num_classes, self.img_size,self.device)
        if split_name == 'test':
            batch_size = len(dataset)
        else:
            batch_size = self.batch_size
        return DataLoader(dataset, batch_size, shuffle=True)
