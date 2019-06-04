from scipy import ndimage as ndi
from skimage.morphology import dilation, ball
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from skimage.measure import regionprops
from copy import deepcopy
import numpy as np



def watershed_segmentation(myMolecule, level, scale_maps, steps):
    imageData = myMolecule.data()
    thresholded = imageData > level
    mask = dilation(thresholded, ball(2))
    labels = watershed(-imageData, connectivity=26, mask=mask)
    local_maxima = peak_local_max(imageData, labels=labels, num_peaks_per_label=1)

    for step in range(steps):
        step_local_maxima = peak_local_max(scale_maps[step].data())
        for point in np.nditer(local_maxima):
            print("Local maxima ", len(local_maxima))
            print(np.linalg.norm(step_local_maxima-point, axis=1).argmin(axis=0))
            local_maxima[i] = np.linalg.norm(step_local_maxima-point, axis=1).argmin(axis=0)
    print(local_maxima)


    
    #Unify regions

    #atoms_distance["distance"] = distance(atomic_structure["position"])
    return labels

def gaussian_step_filtering(myMolecule, step_sigma, steps):
    step_maps = []
    initial_sigma = step_sigma
    for map_id in range(steps):
        newMolecule = deepcopy(myMolecule)
        smoothed = gaussian(newMolecule.data(), sigma=initial_sigma)
        newMolecule.set_data(smoothed)
        step_maps.append(newMolecule)
        initial_sigma+=step_sigma
    return step_maps
