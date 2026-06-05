import os
from glob import glob
import argparse
import gzip
import mrcfile
from metrics import getCorrelation, getRelative_Masks_Overlap
from concurrent.futures import wait
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd
import copy
import json
import metrics
import time
import gc

from skimage.measure import regionprops
from scipy.ndimage import distance_transform_edt, gaussian_filter

from Bio.PDB import PDBParser, PDBIO
import torch
def convert(o):
    if isinstance(o, np.generic): return o.item()  
    raise TypeError


def save_compressed_npy(path, array):
    temp_path = path + '.tmp'
    try:
        with gzip.open(temp_path, 'wb') as f:
            np.save(f, array)
        os.replace(temp_path, path)
    except Exception:
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except Exception:
                pass
        raise


def load_numpy_array(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            return np.load(f, allow_pickle=False)

    arr = np.load(path, allow_pickle=False)
    if isinstance(arr, np.lib.npyio.NpzFile):
        # load first dataset from npz container
        keys = list(arr.files)
        if keys:
            return arr[keys[0]]
    return arr

# Interseccion de mapas simulados de pedazos con original
# Si hay traslape debe anotarse
# Obtiene mapa anotado segun label, tipo float
# Revisa pedazos no asociados, utiliza holgura, hace una pasada
# obtiene stats
# Lo guarda en disco

def annotateSample(map_id, indexes, df, fullness,columns, output_dir):
    map_id = df.at[indexes[0], columns['id']]
    map_path = '/data/ggutierrez/dataset/models/'+map_id+'_resized.mrc'
    annotated_path = os.path.join(output_dir,os.path.basename(map_path).replace('.','_gt.'))
    contourLvl = float(df.at[indexes[0], columns['contourLevel']])
    map_to_annotate =  mrcfile.open(map_path)
    data_map = map_to_annotate.data
    map_mask = data_map >= contourLvl
    result = {}
    result['map_path'] = map_path
    result['contourLevel'] = contourLvl
    result['total'] = np.sum(map_mask)
    # Set to 0 all voxels outside contour level, otherwise fill with a marker
    marker = 10000
    data_map_copy = copy.deepcopy(data_map)
    data_map_copy[np.logical_not(map_mask)] = 0
    data_map_copy[map_mask] = marker
    labels = []
    chain_label_id_dict = {}
    print('Tagging em map {}'.format(os.path.basename(map_path)))
    try:
        for i in indexes:
            segment_path = '/data/ggutierrez/simulated/simulated_chain/sim_'+map_id+'_'+ df.at[i, columns['chain_id']]+'.mrc'
            segment_label = int(float(df.at[i, columns['chain_label']]))
            chain_label_id_dict[df.at[i,columns['chain_label']]] = df.at[i,columns['chain_id']]
            segment_map = mrcfile.open(segment_path)
            segment_mask = segment_map.data >= 0.9
            print("Number of voxels in segment {}".format(np.sum(segment_mask)))
            masks_intersec = np.logical_and(map_mask, segment_mask)
            print("Number of voxels in intersection {}".format(np.sum(masks_intersec)))
            data_map_copy[masks_intersec] = segment_label
            labels.append(segment_label)
            print("Chain {}, voxels {}".format(segment_label,np.sum(segment_mask)))
            print("	Matching {} of {} voxels".format(np.sum(masks_intersec), np.sum(segment_mask)))
            segment_map.close()
    except Exception as e:
            with open('error.txt', 'a') as out:
                out.write('There is a problem getting segments for {}:{}'.format(segment_map,e))
            return ValueError('There is a problem getting segments for {}:{}'.format(segment_map, e))
    #import pdb; pdb.set_trace()
    # Get non assigned voxels
    dim1,dim2,dim3 = np.where(data_map_copy == marker)
    nonassigned_points = np.array(list(map(list,zip(dim1,dim2,dim3))))
    # Get assigned voxels coords
    dim1,dim2,dim3 = np.where(np.logical_and((data_map_copy != marker), (data_map_copy != 0)))
    # Combine list of indexes into a list of points in 3D space
    assigned_points = list(map(list,zip(dim1,dim2,dim3)))
    print("Asigned voxels : {}".format(len(assigned_points)))
    print("Non asigned voxels : {}".format(len(nonassigned_points)))
    print("Total number of voxels: {}".format(np.sum(map_mask)))
    # If any voxel remain
    if (len(nonassigned_points) > 0) & (len(assigned_points)>0):
        print("Attempt to assign {} voxels left".format(len(nonassigned_points)))
        # Create KDTree with assigned points
        tree = cKDTree(assigned_points)
        # Search for nearest point
        d,i = tree.query(nonassigned_points)
        neighbors_index = tree.data[i].astype(int)
        # Use voxels inside fullnes value only
        mask = d <= fullness
        mask_inv = np.logical_not(mask)
        points_to_reassign = nonassigned_points[mask]
        #points_to_reassign = nonassigned_points
        #points_to_discard = []
        points_to_discard = nonassigned_points[mask_inv]
        neighbors_index = neighbors_index[mask]
        d1_i, d2_i, d3_i = neighbors_index[:,0], neighbors_index[:,1], neighbors_index[:,2]
        # Replace values in map with search result
        values_to_map = data_map_copy[d1_i,d2_i,d3_i]
        for point,value in zip(points_to_reassign,values_to_map):
            data_map_copy[point[0],point[1],point[2]] = value
        # Set voxels outside fullness value to 0
        for point in points_to_discard:
            data_map_copy[point[0],point[1],point[2]] = 0
        #pdb.set_trace()
        #result['voxels_reasigned'] = len(points_to_reassign)
        #result['voxels_discarted'] = len(points_to_discard)
        result['voxels_reasigned'] = np.sum(mask)
        result['voxels_discarted'] = np.sum(mask_inv)
    else:
        print("	No more voxels to assign")
        result['voxels_reasigned'] = 0
        result['voxels_discarted'] = 0
    marker_left = np.sum(data_map_copy == marker)
    if marker_left>0:
        print("there shuldnt be {} markers in array of labels.. check this {}".format(marker_left,os.path.basename(map_path)))
    # print labels
    voxels_dict = {}
    for l in labels:
        voxels_dict[l]=np.sum(data_map_copy==l)
        filename = map_path.replace(str(map_path[-4:]), '_'+str(l)+'.npy')
        print("Voxels for label {} :{}".format(l, voxels_dict[l]))
        data_masked = np.copy(data_map_copy)
        data_masked[data_map_copy==l] = 2.0
        data_masked[data_map_copy!=l] = 1.0
        data_masked[data_map_copy==0] = 0.0        
        print("saved volume of {}".format(np.sum(data_masked == 2)))
        np.save(filename, data_masked)
        print("saved {}".format(filename))
        del data_masked
    # Compute euler numbers
    '''
    euler_dict = {}
    for region in regionprops(data_map.astype(np.int32)):
        euler_dict[region.label] = region.euler_number
    # Save map
    result['euler_segments'] = json.dumps(euler_dict, default=convert)
    '''
    result['voxels_assigned'] = json.dumps(voxels_dict, default=convert)
    result['tagged_path'] = annotated_path
    result['map_id'] = map_id
    
    map_gt = mrcfile.open(annotated_path, 'w+')
    map_gt.set_data(data_map_copy)
    
    map_gt.close()
    map_to_annotate.close()
    del data_map_copy

    return result

def generate_sphere_points(n_points=1000):
    print("Generating sphere points...")
    points = []
    phi = np.pi * (3. - np.sqrt(5.))
    
    for i in range(n_points):
        y = 1 - (i / float(n_points - 1)) * 2
        radius = np.sqrt(1 - y * y)
        
        theta = phi * i
        
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        
        points.append([x, y, z])
    
    return np.array(points)

def select_surface_points_from_sphere(region_gt, density_map, contour_level, n_points=1000):
    print("Selecting surface points from sphere...")
    center = np.mean(np.where(region_gt > 0), axis=1)
    
    distance = distance_transform_edt(region_gt)
    distance[distance != 1] = 0
    surface_points = np.array(np.where(distance == 1)).T
    
    # **CRITICAL FIX**: Handle empty surface (e.g., single voxel or no distance==1)
    if len(surface_points) == 0:
        print("  WARNING: No surface voxels found...")
        edge_points = np.array(np.where(region_gt > 0)).T
        if len(edge_points) == 0:
            return np.array([])  # Completely empty region
        surface_points = edge_points  # Fall back to all region points
    
    density_values = density_map[surface_points[:,0], surface_points[:,1], surface_points[:,2]]

    # Handle case where all density values are below contour level
    valid_density = density_values[density_values >= float(contour_level)]
    if len(valid_density) == 0:
        print(f"  WARNING: No density values >= {contour_level}. Using all available densities.")
        percentile_density = np.percentile(density_values, 75) if len(density_values) > 0 else 0
    else:
        percentile_density = np.percentile(valid_density, 75)

    sphere_points = generate_sphere_points(n_points)
    max_radius = np.max(np.linalg.norm(surface_points - center, axis=1))
    sphere_points = sphere_points * max_radius + center
    
    print("Sphere points generated, selecting based on density...")
    selected_points = []
    for sphere_point in sphere_points:
        direction = sphere_point - center
        direction = direction / np.linalg.norm(direction)
        
        distances = np.abs(np.cross(surface_points - center, direction)).sum(axis=1)
        # **CRITICAL FIX**: Check if distances is empty before argmin
        if len(distances) == 0:
            print("  ERROR: Distances array is empty!")
            break
        
        closest_point_idx = np.argmin(distances)
        
        point_density = density_values[closest_point_idx]
        if point_density >= percentile_density:
            selected_points.append(surface_points[closest_point_idx])
    
    print("Selected {} points from sphere based on density.".format(len(selected_points)))
    return np.array(selected_points)

def annotatePoints(df, i, output_path, final_output_path=None, pool_size=1, number_points=3, gaussian_std=1):
    map_path = df.iloc[i]['map_path']
    output_rows = []
    tagged_map_path = df.iloc[i]['tagged_path']
    
    # **FIX**: Ensure mrcfile data is fully loaded into memory with explicit copy
    tagged_map = mrcfile.open(tagged_map_path)
    try:
        # Force full read and copy to memory
        tagged_map_data = np.array(tagged_map.data, dtype=np.float32, copy=True)
    except Exception as e:
        print(f"ERROR: Failed to load mrcfile {tagged_map_path}: {e}")
        tagged_map.close()
        return pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
    
    bbox = None
    distance = None
    target_npy_dir = output_path if output_path is not None else final_output_path
    os.makedirs(target_npy_dir, exist_ok=True)
    
    for region in regionprops(tagged_map_data.astype(np.int32)):
        label = int(region.label)
        bbox = region.bbox
        region_gt = np.copy(tagged_map_data)
        region_gt[ region_gt != label ] = 0.0
        region_gt[ region_gt == label ] = 1.0
        distance = distance_transform_edt(region_gt)
        distance[distance != 1] = 0
        print("creating point set for region {} with volume {} of EM map {}".format(label, np.sum(region_gt), df.iloc[i]['tagged_path'])) 
        if np.sum(region_gt)==0.0:
            print("Tagged path {} does not have assigned voxels, ommiting extreme point anotation".format(df.iloc[i]['tagged_path']))
            continue
        # Fix to include numpy array path
        tagged_path = map_path.replace(str(df.iloc[i]['map_path'][-4:]), '_'+str(label)+'.npz')
        for p in range(pool_size):
            basename = df.iloc[i]['id']+'_' + str(label)+'_' + str(p)+'.npy.gz'
            print("Creating point sample {} for annotated {} ".format(p,basename))
            region_path = os.path.join(target_npy_dir,basename)
            index_x, index_y, index_z = np.where(distance == 1)
            surface_points = select_surface_points_from_sphere(
                region_gt,
                tagged_map_data,
                df.iloc[i]['contourLevel'])
            
            # **FIX**: Validate surface_points before sampling
            if len(surface_points) == 0:
                print(f"WARNING: No surface points found for region {label}. Skipping point sample {p}.")
                continue
            
            if len(surface_points) < number_points:
                print("Warning: only {} surface points available, need {}; using replace=True".format(len(surface_points), number_points))
                chosen_indexes = np.random.choice(len(surface_points), number_points, replace=True)
            else:
                chosen_indexes = np.random.choice(len(surface_points), number_points, replace=False)
            # chosen_indexes = np.random.choice(len(index_x), number_points, replace=False)
            index_x = surface_points[chosen_indexes][:,0]
            index_y = surface_points[chosen_indexes][:,1]
            index_z = surface_points[chosen_indexes][:,2]
            point_array = np.zeros_like(region_gt, dtype=np.float32)
            point_array[index_x,index_y,index_z] = 1.0
            point_array = gaussian_filter(point_array, gaussian_std)
            point_array = point_array.astype(np.float32)
            save_compressed_npy(region_path, point_array)
            output_rows.append({
                'id': df.iloc[i]['id'],
                'map_path': df.iloc[i]['map_path'],
                'contourLevel': df.iloc[i]['contourLevel'],
                'subunit': label,
                'tagged_path': tagged_path,
                'number_points': number_points,
                'tagged_points_path': region_path,
                'min_x': bbox[0],
                'min_y': bbox[1],
                'min_z': bbox[2],
                'max_x': bbox[3],
                'max_y': bbox[4],
                'max_z': bbox[5]
            })
            del point_array
    
    tagged_map.close()
    del tagged_map
    del tagged_map_data
    if distance is not None:
        del distance
    
    if len(output_rows) > 0:
        output_df = pd.DataFrame(output_rows)
    else:
        output_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
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
    sample_map = mrcfile.open(tagged_path)
    sample = sample_map.data >= contour
    labeled_map = mrcfile.open(compare_path)
    labeled = labeled_map.data >= contour
    iou = metrics.intersection_over_union(sample, labeled)
    h = metrics.homogenity(sample, labeled)
    p = metrics.proportion(sample, labeled)
    c = metrics.consistency(sample, labeled)
    sample_map.close()
    labeled_map.close()
    return pd.Series( [map_id, row['map_path'], tagged_path, row['contourLevel'], compare_path, iou, h, p, c ], index=['id', 'map_path','tagged_path', 'contourLevel', 'reference_path', 'iou', 'homogenity', 'proportion', 'consistency'])

def doParallelTagging(df, fullness, gt_path, columns, comm, size ):
    unique_id_list = df[columns['id']].unique().tolist()
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','map_path','contourLevel','tagged_path','subunits','matched_subunits','voxels','voxels_matched','voxels_discarted','voxels_reassigned','voxels_assigned'])
    print("Spawn procecess...")
    '''  
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
                    tagged_path = res['tagged_path']
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
                    output_df = output_df.append({'id':map_id, 'map_path':map_path, 'contourLevel':contour, 'tagged_path':tagged_path, 'subunits':len(voxels_assigned.keys()), 'matched_subunits':segments_matched, 'voxels':voxels_num, 'voxels_matched':voxels_matched, 'voxels_discarted':voxels_discarted, 'voxels_reassigned':voxels_reassigned, 'voxels_assigned':voxels_assigned}, ignore_index=True)
                    
                except ValueError as error:
                    print("Error asignating segments for {}".format(map_id))
    '''
    for i in unique_id_list:
        subunit_indexes = df.loc[df[columns['id']]==i].index.tolist()
        res = annotateSample(i, subunit_indexes, df, fullness, columns, gt_path)
        map_id = res['map_id']
        voxels_assigned = json.loads(res['voxels_assigned'])
        voxels_reassigned = res['voxels_reasigned']
        voxels_discarted = res['voxels_discarted']
        tagged_path = res['tagged_path']
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
        new_row = pd.DataFrame({'id':[map_id], 'map_path':[map_path], 'contourLevel':[contour], 'tagged_path':[tagged_path], 'subunits':[len(voxels_assigned.keys())], 'matched_subunits':[segments_matched], 'voxels':[voxels_num], 'voxels_matched':[voxels_matched], 'voxels_discarted':[voxels_discarted], 'voxels_reassigned':[voxels_reassigned], 'voxels_assigned':[voxels_assigned]})
        output_df = pd.concat([output_df, new_row], ignore_index=True)
     
    return output_df


def samplingPatches(df, size, extra_width, stride):
    id_list = df.index.tolist()
    output_df = pd.DataFrame(columns=['id','subunit','patch','data_path'])
    
    row = df.iloc[0]

    map_id = row['id']
    segment_id = row['subunit']
    map_object = mrcfile.open(row['map_path'])
    map_data = map_object.data
    map_data = np.copy(map_object.data)
    map_data[map_data<float(row['contourLevel'])]=0
    mask_data = load_numpy_array(row['tagged_path'])

    data_max = np.max(map_data)
    data_min = np.min(map_data)
    norm_data = (map_data - data_min)/ (data_max-data_min + 1e-6)
    
    min_x = max(int(row['min_x'])-extra_width, 0 )
    min_y = max(int(row['min_y'])-extra_width, 0)
    min_z = max(int(row['min_z'])-extra_width, 0)
    max_x = min(int(row['max_x'])+extra_width, map_data.shape[0])
    max_y = min(int(row['max_y'])+extra_width, map_data.shape[1])
    max_z = min(int(row['max_z'])+extra_width, map_data.shape[2])

    data_shape = (max_x-min_x,max_y-min_y,max_z-min_z)

    # Slice
    norm_data = norm_data[min_x:max_x,min_y:max_y,min_z:max_z]
    mask_data = mask_data[min_x:max_x,min_y:max_y,min_z:max_z]

    #Check padding
    # Calculate padding to fit the sliding windows
    pad0_left = (data_shape[0] // stride * stride + size) - data_shape[0]
    pad1_left = (data_shape[1] // stride * stride + size) - data_shape[1]
    pad2_left = (data_shape[2] // stride * stride + size) - data_shape[2]
    # Calculate symmetric padding
    pad0_right = pad0_left // 2 if pad0_left % 2 ==0 else pad0_left // 2 + 1
    pad1_right = pad1_left // 2 if pad1_left % 2 ==0 else pad1_left // 2 + 1
    pad2_right = pad2_left // 2 if pad2_left % 2 ==0 else pad2_left // 2 + 1
    pad0_left = pad0_left // 2
    pad1_left = pad1_left // 2
    pad2_left = pad2_left // 2
    print("Map {} chain {}: Sizes {} and {}, Padding:{},{}:{},{}:{}".format(map_id,segment_id,norm_data.shape, pad0_left, pad0_right, pad1_left, pad1_right, pad2_left, pad2_right))
    norm_data_padded = np.pad(norm_data, ((pad0_left, pad0_right), (pad1_left, pad1_right), (pad2_left, pad2_right)))
    mask_data_padded = np.pad(mask_data, ((pad0_left, pad0_right), (pad1_left, pad1_right), (pad2_left, pad2_right)))
    new_shape = norm_data_padded.shape
    all_tensor = np.zeros((12,new_shape[0],new_shape[1],new_shape[2]))
    filename_tmp = 'data_patches/{}_{}.npy'.format(map_id,segment_id)
    all_tensor[0] = mask_data_padded
    all_tensor[1] = norm_data_padded
    df_count = 0
    
    for index,row in df.iterrows():
        point_id =row['tagged_points_path'][-5:-4]
        point_data = load_numpy_array(row['tagged_points_path'])
        point_data = point_data[min_x:max_x,min_y:max_y,min_z:max_z]
        point_data_padded = np.pad(point_data, ((pad0_left, pad0_right), (pad1_left, pad1_right), (pad2_left, pad2_right)))
        all_tensor[df_count+2] = point_data_padded
        
    all_tensor = torch.from_numpy(all_tensor)           
    patches_tensor = all_tensor.unfold(3, size, size//2).unfold(2, size, size//2).unfold(1, size, size//2)
    patches_tensor = patches_tensor.contiguous().view(-1,12,size,size,size) 
    number_patches = patches_tensor.size(0)
    for i in range(number_patches):#['id','subunit','patch','data_path', 'min_x', 'min_y', 'min_z', 'max_x', 'max_y','max_z']
        new_row = pd.DataFrame({'id':[map_id], 'subunit':[segment_id],'patch':[i],'data_path':['data_patches/{}_{}.npy'.format(map_id,segment_id)]})
        output_df = pd.concat([output_df, new_row], ignore_index=True)
    np.save('data_patches/{}_{}.npy'.format(map_id,segment_id), all_tensor.numpy())
    print("Saved id {} subunit {} with shape {}".format(map_id,segment_id, patches_tensor.shape))
    map_object.close()
    del map_data
    del mask_data
    del norm_data
    del point_data
    del all_tensor
    return output_df
                    
        
            
""" def doParallelSamplingPatches(df, shape_size,stride, extra_width, comm, size ):
    unique_df = df.groupby(['id','subunit'])
    # Construct dataframe to store results
    output_df = pd.DataFrame(columns=['id','subunit','patch','data_path'])
    print("Spawn procecess...")
     
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            # For each map, perform annotation
            for name,group in unique_df:
                futures.append(executor.submit(samplingPatches,group, shape_size, extra_width, stride))
            wait(futures)
            for f in futures:
                try:
                    res = f.result()
                    output_df = output_df.append(res, ignore_index=True)
                except ValueError as error:
                    print("Error asignating patches  for {}".format(error))
    ''' 
    for name,group in unique_df:
        out = samplingPatches(group, shape_size, extra_width, stride)
    '''
    return output_df """

def doParallelAdjacency(df):
    id_list = df.index.tolist()
    print("Spawn procecess...")
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
        output_df = pd.concat([output_df, res.to_frame().T], ignore_index=True)
    return output_df

def doParallelExtremePointAnnotation(df, pool_size, output_path, final_path=None):
    import shutil
    indexes = df.index.tolist()
    
    # If final_path not provided, use output_path for final results
    if final_path is None:
        final_path = output_path
    
    # Define checkpoint paths for recovery
    checkpoint_dir = final_path
    aggregate_checkpoint_path = os.path.join(checkpoint_dir, '.aggregate_checkpoints.csv')
    processed_log_path = os.path.join(checkpoint_dir, '.processed_indexes.txt')
    failed_log_path = os.path.join(checkpoint_dir, '.failed_rows.txt')
    final_csv_path = os.path.join(final_path, 'dataset_extreme_points.csv')
    
    def row_key(idx):
        row = df.iloc[idx]
        subunit = row.get('subunit', '') if hasattr(row, 'get') else ''
        return f"{row['id']}|{row['tagged_path']}|{subunit}"

    def write_log(path, line):
        try:
            with open(path, 'a') as f:
                f.write(line + '\n')
                f.flush()
                os.fsync(f.fileno())
        except Exception:
            pass

    def append_checkpoint_rows(rows):
        if rows is None or len(rows) == 0:
            return
        header = not os.path.exists(aggregate_checkpoint_path)
        rows.to_csv(aggregate_checkpoint_path, mode='a', header=header, index=False)
        try:
            with open(aggregate_checkpoint_path, 'a') as f:
                f.flush()
                os.fsync(f.fileno())
        except Exception:
            pass

    # Load previously processed row keys to enable recovery
    processed_keys = set()
    if os.path.exists(processed_log_path):
        try:
            with open(processed_log_path, 'r') as f:
                processed_keys = set(line.strip() for line in f if line.strip())
            print(f"✓ Recovery: Found {len(processed_keys)} previously processed row keys")
        except Exception as e:
            print(f"Warning: Could not load processed log: {e}")
    
    # Load aggregate checkpoint if it exists
    checkpoint_rows = 0
    if os.path.exists(aggregate_checkpoint_path):
        try:
            with open(aggregate_checkpoint_path, 'r') as f:
                header = f.readline()
                checkpoint_rows = sum(1 for _ in f)
            print(f"✓ Recovery: Found {checkpoint_rows} rows in checkpoint")
        except Exception as e:
            print(f"Warning: Could not load checkpoint: {e}")
    
    # Process in batches of 5 indexes (reduced from 10 for lower memory)
    batch_size = 5
    items_to_process = [idx for idx in indexes if row_key(idx) not in processed_keys]
    total_items = len(indexes)
    
    print(f"\nProcessing Summary: Total={total_items}, Done={len(processed_keys)}, Remaining={len(items_to_process)}\n")
    
    for batch_num, batch_start in enumerate(range(0, len(items_to_process), batch_size)):
        batch_end = min(batch_start + batch_size, len(items_to_process))
        batch_indexes = items_to_process[batch_start:batch_end]
        
        # Progress calculation
        current_batch = batch_num + 1
        total_batches = (len(items_to_process) + batch_size - 1) // batch_size
        processed_so_far = len(processed_keys) + batch_start

        print(f"Batch {current_batch}/{total_batches} | Progress: {processed_so_far}/{total_items}")
        print("-" * 100)
        
        # Process each index in the batch
        for i in batch_indexes:
            row_key_value = row_key(i)
            try:
                res = annotatePoints(df, i, output_path, final_output_path=final_path, pool_size=pool_size)

                append_checkpoint_rows(res)
                write_log(processed_log_path, row_key_value)
                processed_keys.add(row_key_value)
                checkpoint_rows += len(res)

                print(f"  ✓ Index {i:4d} processed | Checkpoint: {checkpoint_rows} rows")
                gc.collect()

            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}"
                write_log(failed_log_path, f"{row_key_value}\t{error_msg}")
                print(f"  ✗ Index {i:4d} ERROR: {error_msg[:80]}")
                print(f"     (Saved: {checkpoint_rows} rows in checkpoint)")
                continue
    
    # Consolidate all checkpoints into final dataset
    print("\n" + "=" * 100)
    print("CONSOLIDATING CHECKPOINT INTO FINAL DATASET")
    print("=" * 100)
    
    aggregate_df = None
    if os.path.exists(aggregate_checkpoint_path):
        try:
            final_df = pd.read_csv(aggregate_checkpoint_path, dtype=str)
            if len(final_df) > 0:
                if {'id', 'tagged_path', 'subunit'}.issubset(final_df.columns):
                    final_df['_row_key'] = (
                        final_df['id'].astype(str) + '|' +
                        final_df['tagged_path'].astype(str) + '|' +
                        final_df['subunit'].astype(str)
                    )
                else:
                    final_df['_row_key'] = final_df.index.astype(str)
                final_df = final_df.drop_duplicates(subset=['_row_key']).drop(columns=['_row_key'])
                final_df.to_csv(final_csv_path, index=False)
                aggregate_df = final_df
                print(f"✓ Created final dataset: {len(final_df)} rows")
            else:
                if not os.path.exists(final_csv_path):
                    empty_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
                    empty_df.to_csv(final_csv_path, index=False)
                    aggregate_df = empty_df
                    print(f"✓ Created empty final dataset header: {final_csv_path}")
                else:
                    print(f"✓ Final dataset already exists and checkpoint is empty: {final_csv_path}")
                    aggregate_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
        except Exception as e:
            print(f"⚠ Could not build final dataset from checkpoint: {e}")
            if not os.path.exists(final_csv_path):
                empty_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
                empty_df.to_csv(final_csv_path, index=False)
                aggregate_df = empty_df
                print(f"✓ Created empty final dataset header: {final_csv_path}")
            else:
                aggregate_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
    else:
        empty_df = pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])
        empty_df.to_csv(final_csv_path, index=False)
        aggregate_df = empty_df
        print(f"✓ Created empty final dataset header: {final_csv_path}")

    print(f"\n✓ Final dataset: {final_csv_path}")
    print("=" * 100 + "\n")

    def safe_move_file(src_path, dest_path):
        try:
            dest_dir = os.path.dirname(dest_path)
            if os.path.exists(src_path):
                file_size = os.path.getsize(src_path)
                free_space = shutil.disk_usage(dest_dir).free
                if free_space < file_size + 1024 * 1024:
                    print(f"Skipping move, not enough space for {os.path.basename(src_path)}: need {file_size} bytes, have {free_space} bytes")
                    return False

            if os.path.exists(dest_path):
                os.remove(dest_path)
            shutil.move(src_path, dest_path)
            return True
        except Exception as move_error:
            # Clean up any partial destination file created during move
            try:
                if os.path.exists(dest_path):
                    os.remove(dest_path)
            except Exception:
                pass
            print(f"Could not move {src_path}: {move_error}")
            return False

    # Binary .npy.gz files will remain in the temp directory (output_path)
    # Only the final CSV dataset is saved to final_path
    print(f"\nBinary output files (.npy.gz) kept in temp directory: {output_path}")
    print(f"Final dataset CSV saved to: {final_csv_path}\n")

    # Keep checkpoint files for reference (user can delete if restarting)
    print("Checkpoint files (kept for recovery, delete to restart from scratch):")
    print(f"  - {aggregate_checkpoint_path}")
    print(f"  - {processed_log_path}\n")

    return aggregate_df if aggregate_df is not None else pd.DataFrame(columns=['id','map_path','contourLevel','subunit', 'tagged_path', 'number_points','tagged_points_path','min_x','min_y','min_z','max_x','max_y','max_z'])

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--a', required=False, default=6, help='Annotate data with a fullness parameter')
    parser.add_argument('--s', required=False, default=96, help='Size of the 3d path for corresponding sample') 
    parser.add_argument('--p', required=False, default=1, 
                      help='Number of point sets to generate')
    parser.add_argument('--n', required=False, default=3, 
                      help='Number of points per set')
    parser.add_argument('--result_dir', required=False,
                      help='Output directory for temporary results (e.g., /home)')
    parser.add_argument('--final_dir', required=False,
                      help='Final output directory (e.g., /data). If not provided, uses result_dir')

    opt = parser.parse_args()
    
    try:
        # Set up output directories
        temp_dir = opt.result_dir
        final_dir = opt.final_dir if opt.final_dir else opt.result_dir
        
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        if not os.path.exists(final_dir):
            os.makedirs(final_dir)
        
        print(f"Temporary output directory: {temp_dir}")
        print(f"Final output directory: {final_dir}")
        print(f"Starting processing at {time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Load dataset in chunks to save memory (reduced chunk_size for lower RAM)
        chunk_size = 25
        # Use script directory to find CSV so it works on cluster
        script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_path = os.path.join(script_dir, 'dataset_exp_tagged.csv')

        print(f"Loading dataset from: {csv_path}")

        for chunk_idx, exp_tagged in enumerate(pd.read_csv(csv_path, 
                                                         dtype=str, 
                                                         chunksize=chunk_size)):
            print()
            print(f"Processing chunk {chunk_idx + 1}")
            
            # **FIX**: Reset index for each chunk to prevent out-of-bounds errors in recovery
            # When pd.read_csv loads chunks, indices are 0-based for each chunk
            exp_tagged = exp_tagged.reset_index(drop=True)
            
            # Process chunk with temporary and final directories
            extreme_points_df = doParallelExtremePointAnnotation(
                exp_tagged,
                pool_size=int(opt.p),
                output_path=temp_dir,
                final_path=final_dir
            )
            
            # Clear memory
            del extreme_points_df
            gc.collect()
            
        print(f"Processing completed at {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Final dataset saved to {os.path.join(final_dir, 'dataset_extreme_points.csv')}")
        
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        raise
    
    # comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    
    # results_path = os.path.join(current_dir, opt.result_dir)
    # if not os.path.exists(results_path):
    #     os.makedirs(results_path)
        
    # # Clear any existing intermediate files
    # for f in glob.glob(os.path.join(results_path, 'intermediate_results_*.csv')):
    #     os.remove(f)

    # # Fullness parameter
    # fullness = int(opt.a)
    # # Pool size parameter
    # pool_size= int(opt.p)
    # patch_size= int(opt.s)

    # df_exp_merged = pd.read_csv('dataset_exp_merged.csv', dtype=str)
    # Do parallel computation, one process for each map
    # Get index list to schedule processess 
    # Get id unique values to extract indexes of respective molecule subunits 
    # exp_tagged = doParallelTagging(df_exp_merged, fullness, gt_path, {'id':'id','contourLevel':'contourLevel', 'chain_label':'chain_label','chain_id':'chain_id'}, comm, size)
    # Perform same procedure for simulated data.
    #df_sim = pd.read_csv('dataset_sim_merged.csv')
    #sim_tagged=  doParallelTagging(df_sim, fullness, gt_path, {'id':'entries','map_path':'map_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label','chain_id':'chain_id'})

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

    # exp_tagged = pd.read_csv('dataset_exp_tagged.csv', dtype=str)

    # extreme_points_df = doParallelExtremePointAnnotation(
    #     exp_tagged, 
    #     pool_size=int(opt.p),
    #     output_path=results_path
    # )

    # extreme_points_df.to_csv('dataset_extreme_points.csv', index = False)

    #df = pd.read_csv('dataset_extreme_points.csv', dtype=str)
    #output_df = doParallelSamplingPatches(df, patch_size,patch_size//2, 8, comm, size) 
    #output_df.to_csv('dataset_patches.csv', index=False)
    #sim_tagged = pd.read_csv('dataset_sim_tagged.csv')
    #extreme_points_df = doParallelExtremePointAnnotation(sim_tagged, pool_size, os.path.join(current_dir,'extreme_points/'))
    #extreme_points_df.to_csv('dataset_extreme_points_sim.csv', index = False)
    #selected_df = pd.read_csv('dataset_selected.csv')
    #result_df = doParallelAdjacency(selected_df)
    #result_df.to_csv('dataset_selected_adjacency.csv', index=False)

if __name__ == '__main__':
    main()
