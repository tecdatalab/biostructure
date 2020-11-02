import os
import subprocess
from glob import glob
import argparse
import sys 
from em import molecule
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from concurrent.futures import wait
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd

# Intersección de mapas simulados de pedazos con original
# Si hay traslape debe anotarse
# Obtiene mapa anotado según label, tipo float
# Revisa pedazos no asociados, utiliza holgura, hace una pasada
# obtiene stats
# Lo guarda en disco

def annotateSample(map_id, indexes, df, fullness, output_dir):
    map_path = df.at[indexes[0], 'aligned_path']
    annotated_path = os.path.join(output_dir,map_path.replace('.','_gt.'))
    contourLvl = float(df.at[indexes[0], 'contourLevel'])
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
        segment_path = df.at[i, 'subunit_path']
        if os.path.exists(segment_path):
            segment_label = int(float(df.at[i, 'chain_label']))
            segment_map = molecule.Molecule(segment_path, recommendedContour=0.001)
            segment_mask = segment_map.getContourMasks()[1]
            print("Number of voxels in segment {}".format(np.sum(segment_mask)))
            masks_intersec = np.logical_and(map_mask, segment_mask)
            print("Number of voxels in intersection {}".format(np.sum(masks_intersec)))
            data_map[masks_intersec] = segment_label
            labels.append(segment_label)
            print("Chain {}, voxels {}".format(segment_label,segment_map.getVolume()[1]))
            print("	Matching {} of {} voxels".format(np.sum(masks_intersec), np.sum(segment_mask)))
            print(np.sum(data_map==1))
            print(np.sum(data_map==2))
            print(np.sum(data_map==3))
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
        result['voxels_descarted'] = np.sum(mask_inv)
    else:
        print("	No more voxels to assign")
        result['voxels_reasigned'] = 0
        result['voxels_descarted'] = 0
    # print labels
    voxels_dict = {}
    for l in labels:
        voxels_dict[l]=np.sum(data_map==l)
    # Save map
    result['voxels_assigned'] = str(voxels_dict)	
    result['tag_path'] = annotated_path
    result['map_id'] = map_id
    map_to_annotate.setData(data_map)
    map_to_annotate.save(annotated_path)    
    return result

        
        
        

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
    df_exp = pd.read_csv('dataset_selected.csv')
    #df_sim = pd.read_csv('synthetic_selected.csv')
    # Do parallel computation, one process for each map
    # Get index list to schedule processess 
    # Get id unique values to extract indexes of respective molecule subunits 
    unique_id_list = df_exp.id.unique().tolist()
    # Construct dataframe to store results
    annotation_exp_df = df_exp[['id','fitted_entries', 'chain_id','chain_label']]
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform annotation
            for i in unique_id_list:
                subunit_indexes = df_exp.loc[df_exp['id']==i].index.tolist()
                futures.append(executor.submit(annotateSample,i, subunit_indexes, df_exp, fullness, gt_path))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    map_id = res['map_id']
                    # Do somenthing with res
                    print("Received {}".format(res))
                    #voxels_reassigned = res['voxels_reasigned']
                    #voxels_descarted = res['voxels_descarted']
                    #voxels_assigned = res['voxels_assigned']
                    #annotation_exp_df
                except ValueError as error:
                    print("Error asignating segments for {}".format(map_id))
    '''
    for l in unique_id_list:
        subunit_indexes = df_exp.loc[df_exp['id']==l].index.tolist()
        print(annotateSample(subunit_indexes, df_exp, fullness, gt_path))
        return 0
    ''' 
if __name__ == '__main__':
    main()
