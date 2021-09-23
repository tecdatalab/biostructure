import numpy as np
import torch
from torch.optim import Adam, AdamW, SGD
from torch.utils.data import DataLoader

import ignite.distributed as idist

from SegmentationDataset import SegmentationDataset
from SegmentationUNet import SegmentationUNet
from CustomLoss import CustomLoss

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
        

def set_seed(seed):
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

def get_optimizer(name, parameters, learning_rate, momentum, weight_decay):
    print(name)
    if name=='SGD':
        return SGD(parameters, lr=learning_rate,momentum=momentum, weight_decay=weight_decay)
    elif name=='Adam':
        return Adam(parameters, lr=learning_rate)


class SegmentationAgent:
    def __init__(self, val_percentage, test_ids, num_classes, optimizer, loss_func,
        pool_size, extra_width, batch_size, img_size, csv_path, seed, learning_rate, 
        momentum, weight_decay, depth, device, mixed_precision, idist):
        """
        A helper class to facilitate the training of the model
        """
        self.device = device
        self.seed = seed
        set_seed(self.seed)
        self.randg = torch.Generator() 
        self.num_classes = num_classes
        self.batch_size = batch_size
        self.img_size = img_size
        self.dataframe = pd.read_csv(csv_path)
        self.depth = depth
        self.loss_func = loss_func
        #  Ensure that only local rank 0 split data
        if idist.get_local_rank() > 0:
            idist.barrier()

        train_split, val_split, test_split = self.make_splits(val_percentage, test_ids, seed)
        train_dataset = SegmentationDataset(train_split, self.num_classes, self.img_size, self.randg, extra_width)
        validation_dataset = SegmentationDataset(val_split, self.num_classes, self.img_size, self.randg, extra_width)
        test_dataset = SegmentationDataset(test_split, self.num_classes, self.img_size, self.randg, extra_width)
        if idist.get_local_rank() == 0:
            idist.barrier()

        self.train_loader = idist.auto_dataloader(train_dataset, batch_size=batch_size, shuffle=True, worker_init_fn=np.random.seed(self.seed))
        self.validation_loader = idist.auto_dataloader(train_dataset, batch_size=batch_size, shuffle=False, worker_init_fn=np.random.seed(self.seed))
        self.test_loader = DataLoader(test_dataset, batch_size=2*batch_size, shuffle=False, worker_init_fn=np.random.seed(self.seed), pin_memory=True)

        self.model = idist.auto_model(SegmentationUNet(self.num_classes, self.device, depth=self.depth))
        self.criterion = CustomLoss(loss_func, self.num_classes, None, self.device, mixed_precision)
        self.optimizer = idist.auto_optim(get_optimizer(optimizer, self.model.parameters(), learning_rate, momentum, weight_decay))
        

    def make_splits(self, val_percentage, test_ids, seed=42):
        """
        Split the data into train, validation and test datasets
        :param val_percentage: A decimal number which tells the percentage of
                data to use for validation
        :param test_ids: Id of EM maps to use for testing
        :param seed: Set seed for reproductivity
        :return: tuples of splits
        """
        # Get id list of EM maps in data
        id_list = self.dataframe.groupby('id')['id'].unique().index.tolist()
        # Shuffle list to split data into train, validation and testing
        np.random.shuffle(id_list)
        # Select validation, training and test samples
        val_num = int(val_percentage * len(id_list))
        print("Val EM maps {}".format(val_num))
        test_maps = [test_id for test_id in test_ids if test_id in id_list ]
        samples = [sample_id for sample_id in id_list if sample_id not in test_maps ] 
        validation_maps = samples[:val_num]
        train_maps = samples[val_num:]
        print("Number of EM Maps in total: {}, Training: {}, Validation: {}, Testing: {}".format(len(id_list), len(train_maps), len(validation_maps), len(test_maps)))
        validation_segments = self.dataframe[self.dataframe['id'].isin(validation_maps)]
        train_segments = self.dataframe[self.dataframe['id'].isin(train_maps)]
        test_segments = self.dataframe[self.dataframe['id'].isin(test_maps)].drop_duplicates(subset=['id','subunit'])
        print("Number of Segments in total: Training: {}, Validation: {}, Testing: {}".format(len(train_segments.drop_duplicates(subset=['id','subunit'])), len(validation_segments.drop_duplicates(subset=['id','subunit'])), len(test_segments)))
        print("Training samples: {}".format(train_maps))
        print("Validation samples: {}".format(validation_maps))
        print("Testing samples: {}".format(test_maps))
        return train_segments, validation_segments, test_segments

    def get_dataloader(self, split, pool_size=1, extra_width=0, split_name=None, shuffle=False):
        """
        Create a DataLoader for the given split
        :param split: train split, validation split or test split of the data
        :param pool_size: Number of samples in extreme points pool
        :return: DataLoader
        """
        dataset = SegmentationDataset(split, self.num_classes, self.img_size, self.randg, self.device, extra_width)
        if split_name == 'test':
            batch_size = len(dataset)
        else:
            batch_size = self.batch_size
        return DataLoader(dataset, batch_size, shuffle=shuffle, worker_init_fn=np.random.seed(self.seed))

