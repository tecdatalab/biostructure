from scipy import ndimage as ndi
from skimage.morphology import watershed, local_maxima
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu, gaussian
from skimage.color import label2rgb
from copy import copy
import numpy as np




def watershed_segmentation(myMolecule):
    imageData = myMolecule.data()
    th = threshold_otsu(imageData)
    binary = imageData > th
    distance = ndi.distance_transform_edt(binary)
    local_max = peak_local_max(distance, indices=False, labels=binary, footprint=np.ones((26,26,26)), exclude_border=1)
    markers = ndi.label(local_max, structure=np.ones((3,3,3)))[0] 
    labels = watershed(-distance,markers,mask=binary)
    return labels

def gaussian_smooth(myMolecule, sigma = 1):
    newMolecule = copy(myMolecule)
    smoothed = gaussian(newMolecule.data(), sigma=sigma)
    newMolecule.set_data(smoothed)
    return newMolecule
