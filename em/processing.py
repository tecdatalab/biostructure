from skimage.morphology import dilation, ball
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from skimage.measure import regionprops
from copy import deepcopy
import numpy as np
import math



def watershed_segmentation(myMolecule, level, scale_maps, steps):
    imageData = myMolecule.data()
    threshold = imageData > level
    mask = dilation(threshold, ball(2))
    imageData[imageData<=level]=0
    labels = watershed(-imageData, connectivity=26, mask=mask)
    
    # Space-scale filtering
    local_maxima = peak_local_max(imageData)
    for step in range(steps):
        data = scale_maps[step].data()
        step_local_maxima = peak_local_max(data)
        # Get densities for corresponding local maxima, needed to choose local maxima with steepest ascent density
        #precomputed_data = [ data[tuple(step_local_maxima[i])] for i in range(step_local_maxima.shape[0])]
        for row in range(local_maxima.shape[0]):
            point = local_maxima[row]
            
            id_closest_local_maxima = np.argsort(np.linalg.norm(step_local_maxima-point, axis=1))

            # Order by minimun distance between local maxima at a diferent scale
            # And choose closest local maxima which represent the steepest ascent in density value
            # TODO order by both locality and maximun ascent step in terms of density value
            for local_maxima_id in id_closest_local_maxima:
                voxel_coords = tuple(step_local_maxima[local_maxima_id])
                if data[voxel_coords] > data[tuple(point)]:
                    local_maxima[row] = step_local_maxima[local_maxima_id]
                    break
    # Unify labels
    vals, inverse, count = np.unique(local_maxima, return_inverse=True, return_counts=True, axis=0)
    idx_repeated = np.where(count > 1)[0]
    rows, cols = np.where(inverse == idx_repeated[:, np.newaxis])
    _, inverse_rows = np.unique(rows, return_index=True)
    idx_repeated = np.split(cols, inverse_rows[1:])
    for i in range(len(idx_repeated)):
        new_label = idx_repeated[i].min()+1  #Labels are from 1 to n
        mask = np.isin(labels, idx_repeated[i]+1) ##Labels are from 1 to n
        labels[mask] = new_label
    
    return labels

def gaussian_step_filtering(myMolecule, level, step_sigma, steps):
    step_maps = []
    newMolecule = deepcopy(myMolecule)
    data = newMolecule.data()
    data[data<=level] = 0
    for map_id in range(steps):
        newMolecule = deepcopy(newMolecule) 
        smoothed = gaussian(newMolecule.data(), sigma=step_sigma)
        newMolecule.set_data(smoothed)
        step_maps.append(newMolecule)
        step_sigma +=step_sigma
    return step_maps
