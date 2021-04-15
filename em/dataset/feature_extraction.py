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
import traceback
import random

import json
from json import encoder

from skimage.measure import regionprops
from scipy.ndimage import distance_transform_edt, gaussian_filter

from Bio.PDB import PDBParser, PDBIO

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
    result['map_path'] = map_path
    result['contourLevel'] = contourLvl
    result['total'] = map_to_annotate.getVolume()[1]
    # Set to 0 all voxels outside contour level, otherwise fill with a marker
    marker = 10000
    data_map[np.logical_not(map_mask)] = 0
    data_map[map_mask] = marker
    labels = []
    chain_label_id_dict = {}
    print('Tagging em map {}'.format(os.path.basename(map_path)))
    for i in indexes:
        segment_path = df.at[i, columns['subunit_path']]
        if os.path.exists(segment_path):
            segment_label = int(float(df.at[i, columns['chain_label']]))
            chain_label_id_dict[df.at[i,columns['chain_label']]] = df.at[i,columns['chain_id']]
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
    dim1,dim2,dim3 = np.where(data_map == marker)
    if len(dim1)>0:
        print("there shuldnt be markers in array of labels.. check this {}".format(os.path.basename(map_path)))
    # print labels
    voxels_dict = {}
    for l in labels:
        voxels_dict[l]=np.sum(data_map==l)
        filename = map_path.replace(str(map_path[-4:]), '_'+chain_label_id_dict[l]+'.npy')
        map_masked = np.copy(data_map)
        print("Voxels for label {} :{}".format(l, voxels_dict[l])) 
        map_masked[data_map==l] = 1.0
        map_masked[data_map!=l] = 0.0
        print("saved volume of {}".format(map_masked.sum()))
        np.save(filename, map_masked)
        print("saved {}".format(filename))
    # Compute euler numbers
    euler_dict = {}
    for region in regionprops(data_map.astype(np.int32)):
        euler_dict[region.label] = region.euler_number
    # Save map
    result['euler_segments'] = json.dumps(euler_dict, default=convert)
    result['voxels_assigned'] = json.dumps(voxels_dict, default=convert)
    result['tag_path'] = annotated_path
    result['map_id'] = map_id
    map_to_annotate.setData(data_map)
    map_to_annotate.save(annotated_path)    
    return result

def annotatePoints(df, i, output_path, number_points=3, gaussian_std=3):
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path'])
    #print("aa{}".format(df.iloc[i]['tagged_path']))
    tagged_map = molecule.Molecule(df.iloc[i]['tagged_path'], 0.001).getEmMap().data()
    #print("unique",np.unique(tagged_map))
    for region in regionprops(tagged_map.astype(np.int32)):
        label = int(region.label)
        region_gt = np.copy(tagged_map)
        region_gt[ region_gt != label ] = 0.0
        region_gt[ region_gt == label ] = 1.0
        #print("number",np.sum(region_gt==1.0))
        #print("in label {}".format(label))
        basename = df.iloc[i]['id']+'_'+str(label)+'.npy'
        region_path = os.path.join(output_path,basename)
        #print("pathh {}".format(region_path))
        distance = distance_transform_edt(region_gt)
        distance[distance != 1] = 0
        index_x, index_y, index_z = np.where(distance == 1)
        chosen_indexes = np.random.choice(len(index_x), number_points, replace=False)
        #print("indexes:",chosen_indexes)
        index_x = index_x[chosen_indexes]
        index_y = index_y[chosen_indexes]
        index_z = index_z[chosen_indexes]
        point_array = np.zeros_like(region_gt)
        point_array[index_x,index_y,index_z] = 1.0
        point_array = gaussian_filter(point_array, gaussian_std)
        np.save(region_path,point_array) 
        #print("saved {}".format(np.sum(point_array)))
        output_df = output_df.append({'id':df.iloc[i]['id'], 'map_path':df.iloc[i]['map_path'], 'contourLevel':df.iloc[i]['contourLevel'], 'subunit':label, 'tagged_path':df.iloc[i]['tagged_path'], 'number_points':number_points, 'tagged_points_path':region_path}, ignore_index=True)
    #print("output_df: ", output_df)
    return output_df
        
def compute_adjacency(df, i):
    # Get EM map id
    map_id = df.iloc[i]['id']
    # Get pdb path and chain id 
    pdb_path = df.iloc[i]['pdb_path']
    chain = df.iloc[i]['fitted_entries']
    # Create parser and get readed object
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    pdb_obj = parser.get_structure(chain, pdb_path)
    # Compute dictionary to translate chain id (letter) to chain label (number)
    chain_id_list = [chain._id for chain in pdb_obj.get_chains()]
    chain_label_list = [i for i in range(1,len(chain_id_list)+1)]
    dict_label_id_chain = dict(zip(chain_id_list,chain_label_list))
    # Create dictionaries to store coords and kdtree for each chain
    dict_chain_kdtree = dict()
    # Create dictionary to store final adjency data
    adjacency_dict = dict()
    # Compute kdtree for each chain and assign it along with their coords to the corresponding chain label in dict
    for c in pdb_obj.get_chains():        
        ca_coord_list = [atom.coord for atom in c.get_atoms() if atom.name=="CA"]
        chain_id = c.id
        print("get {} atoms for chain {}".format(len(ca_coord_list), chain_id))
        if len(ca_coord_list) == 0:
            continue
        else:
             kdtree = cKDTree(ca_coord_list)
             dict_chain_kdtree[dict_label_id_chain[chain_id]] = kdtree
    # Loop over chains again to compute adjacency (if exists an atom from other chain at a distance of 4 o less Angstroms )
    for c in dict_chain_kdtree.keys():
        # Get atoms coords for current chain from dict
        current_chain_adjacency_dict = dict()
        current_kdtree = dict_chain_kdtree[c]
        # For every other chain, loop atoms to find adjacency or until atom list is empty.
        for c_i in dict_chain_kdtree.keys():
            if c == c_i:
                continue
            else:
                print("Comparing {} against {}".format(c,c_i))
                # Get kdtree to compare with
                chain_kdtree = dict_chain_kdtree[c_i]
                # Get adjacent atoms within radius of 4 Angstroms
                adjacent_atoms = current_kdtree.query_ball_tree(chain_kdtree, r=5)
                number_adjacencies = np.sum([len(adjacent) for adjacent in adjacent_atoms]) 
                if number_adjacencies > 0:
                    current_chain_adjacency_dict[c_i] = 1
                else:
                    current_chain_adjacency_dict[c_i] = 0
        adjacency_dict[c] = current_chain_adjacency_dict

    label_id_chain = json.dumps(dict_label_id_chain, default=convert)
    adjacency = json.dumps(adjacency_dict, default=convert)

    return pd.Series( [map_id, label_id_chain, adjacency], index=['map_id','chain_id_to_label','adjacency'])          
                    
        
            
    
def mapMetricsCompute(row,match_dict):
    map_id = row['id']
    tagged_path = row['tagged_path']
    contour = 0.001
    compare_path = match_dict[map_id]
    sample = molecule.Molecule(tagged_path, contour)
    labeled = molecule.Molecule(compare_path, contour)
    iou = metrics.intersection_over_union(sample, labeled)
    h = metrics.homogenity(sample, labeled)
    p = metrics.proportion(sample, labeled)
    c = metrics.consistency(sample, labeled)
    return pd.Series( [map_id, row['map_path'], tagged_path, row['contourLevel'], compare_path, iou, h, p, c ], index=['id', 'map_path','tagged_path', 'contourLevel', 'reference_path', 'iou', 'homogenity', 'proportion', 'consistency'])

def doParallelTagging(df, fullness, gt_path, columns):
    unique_id_list = df[columns['id']].unique().tolist()
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','tagged_path','subunits','matched_subunits','voxels','voxels_matched','voxels_discarted','voxels_reassigned','voxels_assigned','euler_segments'])
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
                    euler_segments = json.loads(res['euler_segments'])
                    voxels_reassigned = res['voxels_reasigned']
                    voxels_discarted = res['voxels_discarted']
                    tagged_path = res['tag_path']
                    map_path = res['map_path']
                    contour = res['contourLevel']
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
                    output_df = output_df.append({'id':map_id, 'map_path':map_path, 'contourLevel':contour, 'tagged_path':tagged_path, 'subunits':len(voxels_assigned.keys()), 'matched_subunits':segments_matched, 'voxels':voxels_num, 'voxels_matched':voxels_matched, 'voxels_discarted':voxels_discarted, 'voxels_reassigned':voxels_reassigned, 'voxels_assigned':voxels_assigned, 'euler_segments':euler_segments}, ignore_index=True)
                    
                except ValueError as error:
                    print("Error asignating segments for {}".format(map_id))

    return output_df

def doParallelAdjacency(df):
    id_list = df.index.tolist()
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    output_df = pd.DataFrame(columns=['map_id','chain_id_to_label', 'adjacency'])
    ''' 
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform annotation
            for i in id_list:
                futures.append(executor.submit(compute_adjacency,df,i))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    print("Received {}".format(res))
                    output_df = output_df.append(res, ignore_index=True)
                except Exception as error:
                    print(traceback.format_exc())
    '''
    for i in id_list:
        res = compute_adjacency(df,i)
        output_df = output_df.append(res, ignore_index=True)
    return output_df

def doParallelExtremePointAnnotation(df, output_path):
    indexes = df.index.tolist()
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path'])
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
        # For each map, perform annotation
        for i in indexes:
            futures.append(executor.submit(annotatePoints, df, i, output_path))
        wait(futures)
        for f in futures:
            try:
                res = f.result()
                output_df = output_df.append(res, ignore_index=True)
                print(res)
            except Exception as e:
                print("Error annotating extreme points",e)

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
    #exp_sim_metrics = doParallelTagging(df_exp_merged, fullness, gt_path, {'id':'id','map_path':'simulated_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label', 'chain_id':'chain_id'}) 
    #exp_tagged = doParallelTagging(df_exp_merged, fullness, gt_path, {'id':'id','map_path':'aligned_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label','chain_id':'chain_id'})
    # Perform same procedure for simulated data.
    #df_sim = pd.read_csv('dataset_sim_merged.csv')
    #sim_tagged=  doParallelTagging(df_sim, fullness, gt_path, {'id':'entries','map_path':'simulated_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label','chain_id':'chain_id'})

    #match_exp = dict(zip(exp_tagged['id'], exp_tagged['tagged_path']))
    #match_exp_sim = dict(zip(sim_tagged['id'], sim_tagged['tagged_path']))
    # Compute metrics for each dataframe
    #exp_metrics = exp_tagged.apply(lambda x: mapMetricsCompute(x,match_exp), axis=1)
    #sim_metrics= sim_tagged.apply(lambda x: mapMetricsCompute(x,match_exp_sim), axis=1)

    # Create result dataframe with metrics
    #exp_metrics.to_csv('dataset_exp_metrics.csv', index=False)
    #exp_tagged.to_csv('dataset_exp_tagged.csv', index=False)
    #sim_metrics.to_csv('dataset_sim_metrics.csv', index = False)         
    #sim_tagged.to_csv('dataset_sim_tagged.csv', index=False) 
    
    #extreme_points_df = doParallelExtremePointAnnotation(exp_metrics, os.path.join(current_dir,'extreme_points/'))
    #extreme_points_df.to_csv('dataset_extreme_points.csv', index = False)
    #extreme_points_df = doParallelExtremePointAnnotation(sim_metrics, os.path.join(current_dir,'extreme_points/'))
    #extreme_points_df.to_csv('dataset_extreme_points_sim.csv', index = False)
    selected_df = pd.read_csv('dataset_selected.csv')
    result_df = doParallelAdjacency(selected_df)
    result_df.to_csv('dataset_selected_adjacency.csv', index=False)
if __name__ == '__main__':
    main()
