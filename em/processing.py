from skimage.morphology import dilation, ball
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from copy import deepcopy
import numpy as np



def watershed_segmentation(myMolecule, level, scale_maps, steps):
    imageData = myMolecule.data()
    thresholded = imageData > level
    mask = dilation(thresholded, ball(2))
    labels = watershed(-imageData, connectivity=26, mask=mask)
    local_maxima = peak_local_max(imageData, labels=labels, num_peaks_per_label=1)
    # Space-scale filtering
    for step in range(steps):
        step_local_maxima = peak_local_max(scale_maps[step].data())
        for row in range(local_maxima.shape[0]):
            point = local_maxima[row]
            local_maxima[row] = step_local_maxima[np.linalg.norm(step_local_maxima-point, axis=1).argmin()]
    # Unify labels
    vals, inverse, count = np.unique(local_maxima, return_inverse=True, return_counts=True, axis=0)
    idx_repeated = np.where(count > 1)[0]
    rows, cols = np.where(inverse == idx_repeated[:, np.newaxis])
    _, inverse_rows = np.unique(rows, return_index=True)
    idx_repeated = np.split(cols, inverse_rows[1:])
    print(idx_repeated)
    for i in range(len(idx_repeated)):
        new_label = idx_repeated[i].min()+1  #Labels are from 1 to n
        mask = np.isin(labels, idx_repeated[i]+1) ##Labels are from 1 to n
        labels[mask] = new_label
    return labels

def gaussian_step_filtering(myMolecule, step_sigma, steps):
    step_maps = []
    for map_id in range(steps):
        if map_id == 0:
            newMolecule = deepcopy(myMolecule)
        else:
            newMolecule = deepcopy(newMolecule)
        smoothed = gaussian(newMolecule.data(), sigma=step_sigma)
        newMolecule.set_data(smoothed)
        step_maps.append(newMolecule)
    return step_maps
