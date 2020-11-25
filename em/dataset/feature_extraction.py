import os
import subprocess
from glob import glob
import argparse
import sys 
from em import molecule
from em.dataset import metrics
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from concurrent.futures import wait
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd

import json
from json import encoder


def convert(o):
    if isinstance(o, np.generic): return o.item()  
    raise TypeError

# Intersección de mapas simulados de pedazos con original
# Si hay traslape debe anotarse
# Obtiene mapa anotado según label, tipo float
# Revisa pedazos no asociados, utiliza holgura, hace una pasada
# obtiene stats
# Lo guarda en disco

def annotateSample(map_id, indexes, df, fullness,columns, output_dir):
    map_path = df.at[indexes[0], columns['map_path']]
    annotated_path = os.path.join(output_dir,map_path.replace('.','_gt.'))
    contourLvl = float(df.at[indexes[0], columns['contourLevel']])
    map_to_annotate =  molecule.Molecule(map_path, recommendedContour=contourLvl)
    data_map = map_to_annotate.emMap.data()
    map_mask = map_to_annotate.getContourMasks()[1]
    result = {}
    result['total'] = map_to_annotate.getVolume()[1]
    # Set to 0 all voxels outside contour level, otherwise fill with a marker
    marker = 10000
    data_map[np.logical_not(map_mask)] = 0
    data_map[map_mask] = marker
    labels = []
    print('Tagging em map {}'.format(os.path.basename(map_path)))
    for i in indexes:
        segment_path = df.at[i, columns['subunit_path']]
        if os.path.exists(segment_path):
            segment_label = int(float(df.at[i, columns['chain_label']]))
            segment_map = molecule.Molecule(segment_path, recommendedContour=0.001)
            segment_mask = segment_map.getContourMasks()[1]
            print("Number of voxels in segment {}".format(np.sum(segment_mask)))
            masks_intersec = np.logical_and(map_mask, segment_mask)
            print("Number of voxels in intersection {}".format(np.sum(masks_intersec)))
            data_map[masks_intersec] = segment_label
            labels.append(segment_label)
            print("Chain {}, voxels {}".format(segment_label,segment_map.getVolume()[1]))
            print("	Matching {} of {} voxels".format(np.sum(masks_intersec), np.sum(segment_mask)))
        else:
            return ValueError('There is a problem getting segments for {}'.format(aligned_path))
    # Get non assigned voxels
    dim1,dim2,dim3 = np.where(data_map == marker)
    nonassigned_points = np.array(list(map(list,zip(dim1,dim2,dim3))))
    # Get assigned voxels coords
    dim1,dim2,dim3 = np.where(np.logical_and((data_map != marker), (data_map != 0)))
    # Combine list of indexes into a list of points in 3D space
    assigned_points = list(map(list,zip(dim1,dim2,dim3)))
    print("Asigned voxels : {}".format(len(assigned_points)))
    print("Non asigned voxels : {}".format(len(nonassigned_points)))
    print("Total number of voxels: {}".format(map_to_annotate.getVolume()[1]))
    # If any voxel remain
    if (len(nonassigned_points) > 0) & (len(assigned_points)>0):
        # Create KDTree with assigned points
        tree = cKDTree(assigned_points)
        # Search for nearest point
        d,i = tree.query(nonassigned_points)
        neighbors_index = tree.data[i].astype(int)
        # Use voxels inside fullnes value only
        mask = d <= fullness
        mask_inv = np.logical_not(mask)
        points_to_reassign = nonassigned_points[mask]
        points_to_discard = nonassigned_points[mask_inv]
        neighbors_index = neighbors_index[mask]
        d1_i, d2_i, d3_i = neighbors_index[:,0], neighbors_index[:,1], neighbors_index[:,2]
        # Replace values in map with search result
        values_to_map = data_map[d1_i,d2_i,d3_i]
        for point,value in zip(points_to_reassign,values_to_map):
            data_map[point[0],point[1],point[2]] = value
        # Set voxels outside fullness value to 0
        for point in points_to_discard:
            data_map[point[0],point[1],point[2]] = 0
        result['voxels_reasigned'] = np.sum(mask)
        result['voxels_discarted'] = np.sum(mask_inv)
    else:
        print("	No more voxels to assign")
        result['voxels_reasigned'] = 0
        result['voxels_discarted'] = 0
    # print labels
    voxels_dict = {}
    for l in labels:
        voxels_dict[l]=np.sum(data_map==l)
    # Save map
    result['voxels_assigned'] = json.dumps(voxels_dict, default=convert)
    result['tag_path'] = annotated_path
    result['map_id'] = map_id
    map_to_annotate.setData(data_map)
    map_to_annotate.save(annotated_path)    
    return result

    
def mapMetricsCompute(row,match_dict):
    map_id = row['id']
    tagged_path = row['tagged_path']
    contour = 0
    compare_path = match_dict[map_id]
    sample = molecule.Molecule(s_path, contour)
    labeled = molecule.Molecule(gt_path, contour)
    iou = metrics.intersection_over_union(sample, labeled)
    h = metrics.homogenity(sample, labeled)
    p = metrics.proportion(sample, labeled)
    c = metrics.consistency(sample, labeled)
    row['iou'] = iou
    row['homogenity'] = h
    row['proportion'] = p
    row['consistency'] = c

def doParallelTagging(df, fullness, gt_path, columns):
    unique_id_list = df[columns['id']].unique().tolist()
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','tagged_path','subunits','matched_subunits','voxels','voxels_matched','voxels_discarted','voxels_reassigned','voxels_assigned'])
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform annotation
            for i in unique_id_list:
                subunit_indexes = df.loc[df[columns['id']]==i].index.tolist()
                futures.append(executor.submit(annotateSample,i, subunit_indexes, df, fullness, columns, gt_path))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    map_id = res['map_id']
                    voxels_assigned = json.loads(res['voxels_assigned'])
                    voxels_reassigned = res['voxels_reasigned']
                    voxels_discarted = res['voxels_discarted']
                    tagged_path = res['tag_path']
                    voxels_num = res['total']
                    
                    print("Received {}".format(res))
                    # Get number of segments matched
                    segments_matched = 0
                    voxels_matched = 0
                    for key in voxels_assigned.keys():
                        matched_num = voxels_assigned[key]
                        if matched_num > 0:
                            segments_matched+=1
                            voxels_matched += matched_num
                    #'tagged_path', 'subunits','matched_subunits', 'voxels', 'voxels_matched', 'matched_per_segment'
                    output_df = output_df.append({'id':map_id, 'tagged_path':tagged_path, 'subunits':len(voxels_assigned.keys()), 'matched_subunits':segments_matched, 'voxels':voxels_num, 'voxels_matched':voxels_matched, 'voxels_discarted':voxels_discarted, 'voxels_reassigned':voxels_reassigned, 'voxels_assigned':voxels_assigned}, ignore_index=True)
                    
                except ValueError as error:
                    print("Error asignating segments for {}".format(map_id))

    return output_df

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--a', required=False, default=3, help='Annotate data with a fullness parameter')
    parser.add_argument('--result_dir', default='output', required=False, help='output directory to store results')

    opt = parser.parse_args()
    current_dir = os.getcwd()
    
    results_path = os.path.join(current_dir, opt.result_dir)
    gt_path = os.path.join(current_dir, 'ground_truth')
    if not(os.path.exists(gt_path)):
        os.mkdir(gt_path)
    fullness = int(opt.a)
    df_exp_merged = pd.read_csv('dataset_exp_merged.csv')
    # Do parallel computation, one process for each map
    # Get index list to schedule processess 
    # Get id unique values to extract indexes of respective molecule subunits 
    exp_sim_metrics = doParallelTagging(df_exp_merged, fullness, gt_path, {'id':'id','map_path':'simulated_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label'}) 
    exp_metrics = doParallelTagging(df_exp_merged, fullness, gt_path, {'id':'id','map_path':'aligned_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label'})
    # Perform same procedure for simulated data.
    df_sim = pd.read_csv('dataset_sim_merged.csv')
    sim_metrics=  doParallelTagging(df_sim, fullness, gt_path, {'id':'entries','map_path':'simulated_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label'})

    match_exp = dict(zip(exp_sim_metrics['id'], exp_sim_metrics['tagged_path']))
    match_exp_sim = dict(zip(sim_metrics['id'], sim_metrics['tagged_path']))
    # Compute metrics for each dataframe
    exp_metrics = exp_metrics.apply(mapMetricsCompute, match_exp, axis=1)
    exp_sim_metrics = exp_sim_metrics.apply(mapMetricsCompute, match_exp, axis=1)
    sim_metrics= sim_metrics.apply(mapMetricsCompute, match_exp_sim, axis=1)
    # Create result dataframe with metrics
    exp_metrics = exp_metrics.append(exp_sim_metrics)
    exp_metrics = exp_metrics.append(sim_metrics)
    exp_metrics.to_csv('dataset_metrics.csv')          
    

if __name__ == '__main__':
    main()
