from skimage.morphology import dilation, ball
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from skimage.measure import regionprops
import numpy as np
import math
import time
from functools import partial
from multiprocessing import Pool
import os.path



def segment(emMap, steps, sigma, countourLevels, numContours):
    print(type(countourLevels))
    #voxel_size = tuple(int(i/j) for i,j in zip(emMap.cell_dim(), emMap.grid_size()))
    #parallel segmentation for each contour level
    with Pool(numContours) as p:
        segmentation_func = partial(segmentation_pipeline, emMap=emMap, steps=steps, step_sigma=steps)
        labels_list = p.map(segmentation_func, countourLevels)
    return labels_list   
        

def segmentation_pipeline(contourLevel, emMap, steps, step_sigma):
    image = emMap.data()
    step_maps = scale_space_filtering(image, contourLevel, step_sigma, steps)
    labels = watershed_segmentation(image, contourLevel, step_maps, steps)
    # Not accurate yet
    #data, labels = remove_noise_regions(data, labels, voxel_size, 10000) 
    regions = regionprops(labels)   
    regions = [r for r in regions]
    print("Contour Level: %.4f" % contourLevel)
    print("Number of segmented regions: %d" % len(regions))
    print("    step sigma = %.2f\n    steps = %.2f" % (step_sigma, steps))
    try:
        with open(os.path.join("export/", emMap.name+".txt"), "w") as output_file:
            output_file.write("Number of segmented regions: %d\n" % len(regions))
            output_file.write("    step sigma = %.2f\n    steps = %.2f\n" % (step_sigma, steps))
    except Exception as e:
        print("Could not create output file, ", e)
    return labels


def remove_noise_regions(data, labels, voxel_size, area_th=100):
    t = time.process_time()
    regions = regionprops(labels)
    print(regions[1].image.shape)
    print(np.prod(regions[1].image.shape))
    for r in regions:
        print("region %d" % r.label )
        print("number voxels  %d" % r.area)
        print("shape %s" % (r.image.shape,) )
        print("Volume %d" % np.prod(r.image.shape*voxel_size[0]))
        print("Remove? %s" % ("True" if np.prod(r.image.shape*voxel_size[0])<area_th else "False",) )
    coords_to_remove = [r.coords for r in regions if np.prod(r.image.shape*voxel_size[0])<area_th]
    if len(coords_to_remove)!=0:
        coords_to_remove = np.vstack(coords_to_remove).astype(np.int32)
        data[coords_to_remove[:,0],coords_to_remove[:,1],coords_to_remove[:,2]]=0
        labels[coords_to_remove[:,0],coords_to_remove[:,1],coords_to_remove[:,2]]=0
    return data, labels
    

    

def watershed_segmentation(data, level, scale_maps, steps):
    t = time.process_time()
    threshold = data >= level
    mask = dilation(threshold)
    #data[data<level]=0
    preprocess_time = time.process_time() - t
    t = time.process_time()
    labels = watershed(-data, connectivity=26, mask=mask)
    watershed_time = time.process_time() - t
    regions = regionprops(labels)
    # Space-scale filtering
    t = time.process_time()
    local_maxima = peak_local_max(data)  
    initial_maxima_time = time.process_time() - t

    t = time.process_time()
    region_maxima = {}

    #match maxima to each region
    for point in local_maxima:
        region_maxima[labels[tuple(point)]]=point

    #order maxima list
    ordered_keys = np.sort([*region_maxima])

    local_maxima = []
    # Sort local maxima array to match region order
    for key,i in zip(ordered_keys,range(len(regions))):
        #print(key," ",i)
        local_maxima.append(region_maxima[key])

    region_sort_time = time.process_time() - t

    # Visit neighbor voxel at each dimension
    n_range =[-1,0,1]    
    size = data.shape
    t = time.process_time()
    for step in range(steps):
        #print("step ",step)
        smoothed_data = scale_maps[step]

        new_maxima = np.empty_like(local_maxima)

        for n in range(len(local_maxima)):
            
            start_point = local_maxima[n]
            max_val = smoothed_data[tuple(start_point)]
            point = start_point

            # Greedy search
            while(True):
                for k in n_range:
                    k_index = point[0] + k
                    
                    # Explore only valid voxels in z dimension
                    if k_index >= 0 and k_index<size[0]:
                        for j in n_range:
                            j_index = point[1] + j
                            # Explore only valid voxels in y dimension
                            if j_index >= 0 and j_index<size[1]:
                                for i in n_range:
                                    i_index = point[2] + i
                                    # Explore only valid voxels in x dimension
                                    if i_index >= 0 and i_index<size[2]:
                                        # Get smoothed density value at given coordinate
                                        val = smoothed_data[tuple([k_index,j_index,i_index])]
                                        # if density value is bigger means new local maxima is found
                                        if val > max_val:
                                            max_val = val
                                            point = np.array([k_index, j_index, i_index])
                # If new point is equal to the start one, its over. 
                if np.all(point == start_point):
                    break
                # 
                else:
                    start_point = point

            new_maxima[n] = start_point
        #unify labels here, average of each grouped region
        # Unify labels
        vals, inverse, count = np.unique(new_maxima, return_inverse=True, return_counts=True, axis=0)
        idx_repeated = np.where(count > 1)[0]
        rows, cols = np.where(inverse == idx_repeated[:, np.newaxis])
        _, inverse_rows = np.unique(rows, return_index=True)
        idx_repeated = np.split(cols, inverse_rows[1:])

        for group in idx_repeated:
            for id in group:
                local_maxima[id] = new_maxima[id]
        
        
        

    scale_space_time = time.process_time() - t

    t = time.process_time()

    max_label = len(regions)

    vals, inverse, count = np.unique(local_maxima, return_inverse=True, return_counts=True, axis=0)
    idx_repeated = np.where(count > 1)[0]
    rows, cols = np.where(inverse == idx_repeated[:, np.newaxis])
    _, inverse_rows = np.unique(rows, return_index=True)
    idx_repeated = np.split(cols, inverse_rows[1:])

    repeated_maxima_time = time.process_time() - t


    t = time.process_time()
    '''
    replace_dict = {l.label:l.label for l in regions}

    for label in idx_repeated:
        max_label = max_label+1  
        for l in label:
            replace_dict[l] = max_label

    rep_keys, rep_vals = np.array(list(zip(*sorted(replace_dict.items()))))
    indices = np.digitize(labels, rep_keys, right=True)
    labels = rep_vals[indices]

    newlabels = regionprops(labels)

    '''
    for label in idx_repeated:
        max_label = max_label+1  
        for l in label:
            mask = labels == l+1
            labels[mask]=max_label
    newlabels = regionprops(labels)
    for i,l in enumerate(newlabels):
        mask = labels == l.label
        labels[mask]=i+1
        
    for i,l in enumerate(newlabels):
        mask = labels == l.label
        labels[mask]=i+1
    
    
    label_rename_time = time.process_time() - t

    print(" Preprocessing time ", preprocess_time)
    print(" Watershed time ", watershed_time)
    print(" Initial maxima computation time ",initial_maxima_time)
    print(" Region sorting time ", region_sort_time)
    print(" Scale space time ", scale_space_time)
    print(" Repeated maxima search time ", repeated_maxima_time)
    print(" Labels renaming time ", label_rename_time)

    
    return labels



def scale_space_filtering(data, level, step_sigma, steps):
    scaled_data = []
    for map_id in range(steps):
        smoothed = gaussian(data, sigma=[step_sigma,step_sigma,step_sigma], truncate=1)
        scaled_data.append(smoothed)
        step_sigma +=step_sigma
    return scaled_data


