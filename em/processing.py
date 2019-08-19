from skimage.morphology import dilation, ball
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from skimage.measure import regionprops
from copy import deepcopy
import numpy as np
import math



def watershed_segmentation(myMolecule, level, scale_maps, steps):
    data = myMolecule.data()
    threshold = data > level
    mask = dilation(threshold)
    data[data<=level]=0
    labels = watershed(-data, connectivity=26, mask=mask)
    regions = regionprops(labels)
    # Space-scale filtering
    local_maxima = peak_local_max(data)  

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

    # Visit neighbor voxel at each dimension
    n_range =[-1,0,1]    
    size = myMolecule.shape()
    
    for step in range(steps):
        #print("step ",step)
        smoothed_data = scale_maps[step].data()

        new_maxima = np.empty_like(local_maxima)

        for n in range(len(regions)):
            
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
        
        
        



    max_label = len(regions)

    vals, inverse, count = np.unique(local_maxima, return_inverse=True, return_counts=True, axis=0)
    idx_repeated = np.where(count > 1)[0]
    rows, cols = np.where(inverse == idx_repeated[:, np.newaxis])
    _, inverse_rows = np.unique(rows, return_index=True)
    idx_repeated = np.split(cols, inverse_rows[1:])

    #print(idx_repeated)

    for label in idx_repeated:
        max_label = max_label+1  
        for l in label:
            labels[labels==l+1]=max_label
    newlabels = regionprops(labels)
    for i,l in enumerate(newlabels):
        labels[labels==l.label]=i+1
    
    return labels

    
    

def scale_space_filtering(myMolecule, level, step_sigma, steps):
    step_maps = []
    for map_id in range(steps):
        newMolecule = deepcopy(myMolecule) 
        smoothed = gaussian(myMolecule.data(), sigma=[step_sigma,step_sigma,step_sigma], truncate=1)
        newMolecule.set_data(smoothed)
        step_maps.append(newMolecule)
        step_sigma +=step_sigma
    return step_maps
