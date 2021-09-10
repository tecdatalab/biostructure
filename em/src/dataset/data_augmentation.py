from em import molecule
from em.dataset import metrics
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from concurrent.futures import wait
import os
import argparse
import numpy as np
import pandas as pd
import copy

import json
from json import encoder

from skimage.measure import regionprops
import traceback

def convert(o):
    if isinstance(o, np.generic): return o.item()  
    raise TypeError

def overlapAndSplit(i, steps_range, sigma_range, resample_n, df, offset):
    gt_entry = df.iloc[i]
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','gt_path','gt_subunits','matched_subunits','step','sigma','voxels','voxels_assigned','euler_segments','iou', 'proportion', 'consistency', 'homogenity','label'])
    try:
        molecule_object = molecule.Molecule(gt_entry['map_path'], gt_entry['contourLevel'])
        gt_object = molecule.Molecule(gt_entry['tagged_path'], 0.001)
        molecule_density = molecule_object.getDataAtContour(1)
        gt_labels = gt_object.getDataAtContour(1)
        number_gt_segments = gt_entry['matched_subunits']
        # Remove noise
        molecule_density[gt_labels==0] = 0
        molecule_object.setData(molecule_density)
        # First append ground truth to result dataframe
        iou_gt = metrics.intersection_over_union(gt_object,gt_object)
        proportion_gt = metrics.proportion(gt_object, gt_object)
        homogenity_gt = metrics.homogenity(gt_object, gt_object)
        consistency_gt = metrics.consistency(gt_object, gt_object)

        output_df = output_df.append({'id':gt_entry['id'],'map_path':gt_entry['map_path'], 'contourLevel':gt_entry['contourLevel'], 'gt_path':gt_entry['tagged_path'], \
                                  'gt_subunits':gt_entry['subunits'], 'matched_subunits':gt_entry['matched_subunits'], 'step':0, 'sigma':0, 'voxels':gt_entry['voxels'],\
                                  'voxels_assigned':gt_entry['voxels_assigned'], 'euler_segments':gt_entry['euler_segments'], 'iou':iou_gt, 'proportion':proportion_gt, \
                                  'consistency':consistency_gt, 'homogenity':homogenity_gt,'label':'good'},ignore_index=True)

        # Compute diferent segmentation results with watershed
        for i in range(resample_n):
            segmented_count = 1 
            it = 1
            while(segmented_count==1):
                step = np.random.choice(range(steps_range[0],steps_range[1]+1),1, p=[0.4,0.4,0.2])
                sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
                print("iteration {} on molecule {} with {} steps and {} sigma".format(it,gt_entry['id'], step[0], sigma[0]))
                molecule_object.generateSegments(step[0], sigma[0])
                labels = molecule_object.labels.astype(np.int32)
                label_props = regionprops(labels)
                segmented_count = len(label_props)
                it+=1
            # Get gt labels and random choose one to split
            gt_label_props = regionprops(gt_labels.astype(np.int32))
            gt_label_list = [ l.label for l in gt_label_props ]
            label_can_be_splitted = False
            count = 0
            while(label_can_be_splitted==False):
                label_to_be_splitted = np.random.choice(gt_label_list)
                label_mask = (gt_labels == label_to_be_splitted)
                labels_found = np.unique(labels[label_mask])
                number_segments = len(labels_found)
                if ((number_segments > 1) & (number_segments<60)):
                    print("label {} can be splitted in {} segments for molecule {} sample {} after {} iterations".format(label_to_be_splitted,number_segments,gt_entry['id'],i,it))
                    label_can_be_splitted = True
                if count > len(gt_label_list):
                    step = np.random.choice(range(steps_range[0],steps_range[1]+1),1, p=[0.4,0.4,0.2])
                    sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
                    print("Recomputing iteration {} on molecule {} with {} steps and {} sigma".format(it,gt_entry['id'], step[0], sigma[0]))
                    molecule_object.generateSegments(step[0], sigma[0])
                    labels = molecule_object.labels.astype(np.int32)
                    count = 0
                count += 1
            print("spliting label in {} segments with labels {}".format(number_segments, labels_found))
            np.random.shuffle(labels_found)
            rename_label_dict = {}
            count = len(gt_label_list)
            for l in labels_found:
                if count==len(gt_label_list):
                    rename_label_dict[l]=label_to_be_splitted
                    count+=1
                else:
                    rename_label_dict[l] = count
                    count+=1
            print("Rename label dict {}".format(rename_label_dict))
       
            new_labels_object = copy.deepcopy(gt_object)
            new_labels = copy.deepcopy(gt_labels)
            # Split and assign
            for key in np.sort(list(rename_label_dict.keys())):
                mask = np.logical_and(labels==key, new_labels==label_to_be_splitted)
                print("Assigning label {} to {} voxels from gt".format(rename_label_dict[key],np.sum(mask)))
                new_labels[mask] = rename_label_dict[key]
            new_labels_object.setData(new_labels)
            segment_voxels_dict = {}
            segment_euler_dict = {}
            iou = metrics.intersection_over_union(new_labels_object, gt_object)
            proportion = metrics.proportion(new_labels_object, gt_object)
            consistency = metrics.consistency(new_labels_object, gt_object)
            homogenity = metrics.homogenity(new_labels_object, gt_object)
            splitted_labels_props = regionprops(new_labels.astype(np.int32))
            for l in splitted_labels_props:
                segment_voxels_dict[l.label] = np.sum(new_labels == l.label)
                segment_euler_dict[l.label] = l.euler_number
            dict_to_append = {'id':gt_entry['id'], 'map_path':gt_entry['map_path'], 'contourLevel':gt_entry['contourLevel'], 'gt_path':gt_entry['tagged_path'], 'gt_subunits':gt_entry['subunits'], \
            'matched_subunits':len(splitted_labels_props),'step':int(round(step[0])), 'sigma':sigma[0], 'voxels':gt_entry['voxels'], 'voxels_assigned':json.dumps(segment_voxels_dict,default=convert), 'euler_segments':json.dumps(segment_euler_dict, default=convert),'iou':iou, 'proportion':proportion,\
            'consistency':consistency, 'homogenity':homogenity,'label':'good'}
            print(dict_to_append)
            output_df = output_df.append(dict_to_append, ignore_index=True)
    except Exception as e:
        print("Error computing good segmentation for {}: {}".format(gt_entry['id'],e))
        print(traceback.format_exc())
    return output_df

def applyWatershed(i, steps_range, sigma_range, resample_n, df):
    gt_entry = df.iloc[i]
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','gt_path','gt_subunits','matched_subunits','step','sigma','voxels','voxels_assigned','euler_segments','iou', 'proportion', 'consistency', 'homogenity','label'])
 
    try:
        molecule_object = molecule.Molecule(gt_entry['map_path'], gt_entry['contourLevel'])
        gt_object = molecule.Molecule(gt_entry['tagged_path'], 0.001)
        molecule_density = molecule_object.getDataAtContour(1)
        gt_labels = gt_object.getDataAtContour(1)
        # Remove noise
        molecule_density[gt_labels==0] = 0
        molecule_object.setData(molecule_density)
        # Compute diferent segmentation results with watershed
        for i in range(resample_n):
            segmented_count = 61
            while((segmented_count>60) | (segmented_count==1)):
                step = np.random.choice(range(steps_range[0],steps_range[1]+1),1)
                sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
                molecule_object.generateSegments(step[0], sigma[0])
                labels = molecule_object.labels.astype(np.int32)
                label_props = regionprops(labels)
                segmented_count = len(label_props)
            new_labels_object = copy.deepcopy(gt_object)
            new_labels_object.setData(molecule_object.labels)
            segment_voxels_dict = {}
            segment_euler_dict = {}
            iou = metrics.intersection_over_union(new_labels_object, gt_object)
            proportion = metrics.proportion(new_labels_object, gt_object)
            consistency = metrics.consistency(new_labels_object, gt_object)
            homogenity = metrics.homogenity(new_labels_object, gt_object)
            for l in label_props:
                segment_voxels_dict[l.label] = np.sum(labels == l.label)
                segment_euler_dict[l.label] = l.euler_number
            dict_to_append = {'id':gt_entry['id'], 'map_path':gt_entry['map_path'], 'contourLevel':gt_entry['contourLevel'], 'gt_path':gt_entry['tagged_path'], 'gt_subunits':gt_entry['subunits'], \
                'matched_subunits':segmented_count,'step':int(round(step[0])), 'sigma':sigma[0], 'voxels':gt_entry['voxels'], 'voxels_assigned':json.dumps(segment_voxels_dict,default=convert), 'euler_segments':json.dumps(segment_euler_dict, default=convert),'iou':iou, 'proportion':proportion,\
                'consistency':consistency, 'homogenity':homogenity,'label':'bad'}
            print(dict_to_append)
            output_df = output_df.append(dict_to_append, ignore_index=True)
    except Exception as e:
        print("Error computing bad segmentation for {}: {}".format(gt_entry['id'],e))
        print(traceback.format_exc())
    return output_df



def generateBadSegmentation(df, steps_range, sigma_range, resample_n):
    id_list = df.index.tolist()
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','gt_path','gt_subunits','matched_subunits','step', 'sigma', 'voxels','voxels_assigned','euler_segments','iou', 'proportion', 'consistency', 'homogenity', 'label'])
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform operation
            for i in id_list:
                futures.append(executor.submit(applyWatershed, i, steps_range, sigma_range, resample_n, df))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    
                    print("Received {}".format(res))

                    output_df = output_df.append(res, ignore_index=True)
  
                except ValueError as error:
                    print("Error computing bad segments")
    
    return output_df

def generateGoodSegmentation(df, steps_range, sigma_range, resample_n):
    id_list = df.index.tolist()
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','gt_path','gt_subunits','matched_subunits','step', 'sigma', 'voxels','voxels_assigned','euler_segments','iou', 'proportion', 'consistency', 'homogenity','label'])
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform operation
            for i in id_list:
                futures.append(executor.submit(overlapAndSplit, i, steps_range, sigma_range, resample_n, df, 4))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()

                    print("Received {}".format(res))

                    output_df = output_df.append(res, ignore_index=True)

                except ValueError as error:
                    print("Error computing bad segments")
    return output_df

def main():


    tagged_df = pd.read_csv('dataset_exp_tagged.csv')
    # Only entries with matched subunits
    groundt_df = tagged_df[tagged_df['matched_subunits']>0]
    groundt_df = groundt_df.reset_index()

    #segmented_df = generateBadSegmentation(groundt_df, [2,4], [1,3], 10)
    #segmented_df.to_csv('segmented_bad.csv', index = False) 
    segmented_df = generateGoodSegmentation(groundt_df, [4,6],[1,4],10)
    segmented_df.to_csv('segmented_good.csv', index = False) 

if __name__ == '__main__':
    main()
