from scipy import ndimage as ndi
from skimage.morphology import watershed, erosion, dilation
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu, gaussian
from skimage.color import label2rgb
from copy import copy
import numpy as np

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt



def watershed_segmentation(myMolecule):
    imageData = myMolecule.data()
    th = threshold_otsu(imageData)
    thresholded = imageData > th
    kernel = np.ones((3,3,3))
    eroded = erosion(thresholded, kernel)
    dilated = dilation(eroded, kernel)   
    distance = ndi.distance_transform_edt(dilated)
    distance =  distance/np.linalg.norm(distance.ravel(), np.inf, keepdims=True)
    markers,num = ndi.label(np.where(distance>0.7, 1.0, 0.0), structure=kernel)
    labels = watershed(-distance, markers)
    return labels

def gaussian_smooth(myMolecule, sigma = 1):
    newMolecule = copy(myMolecule)
    smoothed = gaussian(newMolecule.data(), sigma=sigma)
    newMolecule.set_data(smoothed)
    return newMolecule
