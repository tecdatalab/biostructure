import numpy as np
import torch
from torch.optim import Adam, AdamW, SGD
from torch.utils.data import DataLoader

import ignite.distributed as idist

from SegmentationDataset import SegmentationDataset
from SegmentationUNet import SegmentationUNet
from CustomLoss import CustomLoss

from sklearn.model_selection import KFold, train_test_split

import pandas as pd

def set_seed(seed):
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)

def get_optimizer(name, parameters, learning_rate, momentum, weight_decay):
    if name=='SGD':
        return SGD(parameters, lr=learning_rate,momentum=momentum, weight_decay=weight_decay)
    elif name=='Adam':
        return Adam(parameters, lr=learning_rate)


class SegmentationAgent:
    def __init__(self, test_ids, num_classes, num_folds, current_fold, optimizer, loss_func,
        extra_width, batch_size, img_size, csv_path, seed, learning_rate, 
        momentum, weight_decay, depth, device, mixed_precision, alpha, beta, idist):
        '''
        A helper class to facilitate the training of the model
        '''

        self.device = device
        self.seed = seed
        set_seed(self.seed)
        self.randg = torch.Generator() 
        self.randg.manual_seed(self.seed)
        self.num_classes = num_classes
        self.num_folds = num_folds
        self.current_fold = current_fold
        self.batch_size = batch_size
        self.img_size = img_size
        self.dataframe = pd.read_csv(csv_path)
        self.depth = depth
        self.loss_func = loss_func
        self.alpha = alpha
        self.beta = beta
        #  Ensure that only local rank 0 split data
        if idist.get_local_rank() > 0:
            idist.barrier()

        train_split, val_split, test_split = self.make_splits(num_folds, test_ids, seed)
        train_dataset = SegmentationDataset(train_split, self.num_classes, self.img_size, self.randg, self.device, augmentate=True, extra_width=extra_width)
        validation_dataset = SegmentationDataset(val_split, self.num_classes, self.img_size, self.randg, self.device, augmentate=False, extra_width=extra_width, is_validation=True)
        test_dataset = SegmentationDataset(test_split, self.num_classes, self.img_size, self.randg, self.device, augmentate=False, extra_width=extra_width, is_validation=True)
        if idist.get_local_rank() == 0:
            idist.barrier()

        self.train_loader = idist.auto_dataloader(train_dataset, batch_size=self.batch_size, shuffle=True, worker_init_fn=np.random.seed(self.seed))
        self.validation_loader = idist.auto_dataloader(train_dataset, batch_size=self.batch_size, shuffle=False, worker_init_fn=np.random.seed(self.seed))
        self.test_loader = DataLoader(test_dataset, batch_size=self.batch_size, shuffle=False, worker_init_fn=np.random.seed(self.seed), pin_memory=True)

        self.model = idist.auto_model(SegmentationUNet(self.num_classes, self.device, in_channels=2, depth=self.depth))
        self.criterion = CustomLoss(loss_func, self.num_classes, self.device, self.alpha, self.beta,  mixed_precision)
        self.optimizer = idist.auto_optim(get_optimizer(optimizer, self.model.parameters(), learning_rate, momentum, weight_decay))
        

    def make_splits(self, num_folds, test_ids, nfolds, seed=42):
        """
        Split the data into train, validation and test datasets
        :param val_percentage: A decimal number which tells the percentage of
                data to use for validation
        :param test_ids: Id of EM maps to use for testing
        :param seed: Set seed for reproductivity
        :return: tuples of splits
        """
        # Get id list of EM maps in data
        df_subunits = self.dataframe.drop_duplicates(subset=['id','subunit'])
        id_list = df_subunits.index.tolist()
        # Generate folds for cross validation
        splits = KFold(n_splits=num_folds,shuffle=True,random_state=seed)
        # Select testing maps from dataset
        if len(test_ids)>0:
            test_maps = [test_id for test_id in test_ids if test_id in id_list ]
            samples = [sample_id for sample_id in id_list if sample_id not in test_maps ]
        else:
            samples, test_maps = train_test_split(id_list, test_size=0.2, random_state=seed) 
            
        # Generate folds for cross validation
        splits = KFold(n_splits=num_folds,shuffle=True,random_state=seed)
        train_idx = []
        val_idx = []
        for t,v in splits.split(samples):
            train_idx.append(t)
            val_idx.append(v)
        validation_maps = [ samples[i] for i in val_idx[self.current_fold] ] 
        train_maps = [ samples[i] for i in train_idx[self.current_fold] ]
        print("Number of samples in total: {}, Training: {}, Validation: {}, Testing: {}".format(len(id_list), len(train_maps), len(validation_maps), len(test_maps)))
        validation_segments = self.dataframe[self.dataframe['id'].isin(validation_maps)]
        train_segments = self.dataframe[self.dataframe['id'].isin(train_maps)]
        test_segments = self.dataframe[self.dataframe['id'].isin(test_maps)].drop_duplicates(subset=['id','subunit'])
        #print("Number of Segments in total: Training: {}, Validation: {}, Testing: {}".format(len(train_segments.drop_duplicates(subset=['id','subunit'])), len(validation_segments.drop_duplicates(subset=['id','subunit'])), len(test_segments)))
        print("Training samples: {}".format(train_maps))
        print("Validation samples: {}".format(validation_maps))
        print("Testing samples: {}".format(test_maps))
        return train_segments, validation_segments, test_segments


