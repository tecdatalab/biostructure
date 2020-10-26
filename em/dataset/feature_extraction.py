import os
import subprocess
from glob import glob
import argparse
import sys 
from em import molecule
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from scipy.spatial import cKDTree
import numpy as np

# Intersección de mapas simulados de pedazos con original
# Si hay traslape debe anotarse
# Obtiene mapa anotado según label, tipo float
# Revisa pedazos no asociados, utiliza holgura, hace una pasada
# obtiene stats
# Lo guarda en disco

def annotateSample(indexes, df, fullness, output_dir):
    map_path = df.at[indexes[0], 'map_path']
    contourLvl = float(df.at[indexes[0], 'contourLevel'])
    map_to_annotate =  molecule.Molecule(map_path, recommendedContour=contourLvl)
    data_map = map_to_annotate.emMap.data()
    map_mask = map_to_annotate.getContourMasks()[1]
    # Set to 0 all voxels outside contour level, otherwise fill with a marker
    marker = 10000
    data_map[np.logical_not(map_mask)] = 0
    data_map[map_mask] = marker
    result = {}
    result['total'] = map_to_annotate.getVolume()[1] 
    labels = []
    for i in indexes:
        print(i)
        segment_path = df.at[i, 'subunit_path']
        segment_label = int(float(df.at[i, 'chain_label']))
        segment_map = molecule.Molecule(segment_path, recommendedContour=1)
        segment_mask = segment_map.getContourMasks()[1]
        masks_intersec = np.logical_and(map_mask, segment_mask)
        data_map[masks_intersec] = segment_label
        labels.append(segment_label)
    # Get non assigned voxels
    dim1,dim2,dim3 = np.where(data_map == marker)
    nonasigned_points = list(map(list,zip(dim1,dim2,dim3)))
    if len(nonasigned_points) > 0 :
        # Get assigned voxels coords
        dim1,dim2,dim3 = np.where((data_map != marker) & (data_map != 0))
        # Combine list of indexes into a list of points in 3D space
        asigned_points = list(map(list,zip(dim1,dim2,dim3)))
        # Create KDTree with assigned points 
        tree = cKDTree(asigned_points)
        # Search for nearest point
        _,i = tree.query(nonasigned_points)
        neighbors_index = np.array(i)
        d1_i, d2_i, d3_i = neighbors_index[:,0], neighbors_index[:,1], neighbors_index[:,2]
        # Replace values in map with search result
        values_to_map = map_data[d1_i][d2_i][d3_i]
        data_map[asigned_points] = values_to_map
    # print labels
    for l in labels:
        result[l]=np.sum(data_map==l)
    # Save map
    
    return result

        
        
        

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--a', required=False, default=3, help='Annotate data with a fullness parameter')
    parser.add_argument('--result_dir', default='output', required=False, help='output directory to store results')

    opt = parser.parse_args()
    current_dir = os.getcwd()
    
    results_path = os.path.join(current_dir, opt.result_dir)
    gt_path = os.path.join(current_dir, 'groundtruth')
    fullness = int(opt.a)
    df_exp = pd.read_csv('dataset_selected.csv')
    #df_sim = pd.read_csv('synthetic_selected.csv')
    # Do parallel computation, one process for each map
    # Get index list to schedule processess 
    # Get id unique values to extract indexes of respective molecule subunits 
    unique_id_list = df_exp.id.unique().tolist()
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform annotation
            for i in uniqueid_list:
                subunit_indexes = df_exp.loc[df_exp['id']==i].index.tolist()
                futures.append(executor.submit(annotateSample, subunit_indexes, fullness, df_exp, gt_path))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    # Do somenthing with res
                except ValueError as error:
                    print("Error calculating volume for simulated {}: {}".format(df_volume.loc[i,'fitted_entries'],error))


if __name__ == '__main__':
    main()
