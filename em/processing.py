from scipy import ndimage as ndi
from skimage.morphology import watershed, erosion, dilation, ball
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from copy import copy
import numpy as np



def watershed_segmentation(myMolecule, level):
    imageData = myMolecule.data()
    thresholded = imageData > level
    mask = dilation(thresholded)

    
    #-------
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    #------
    
    labels = watershed(-imageData, connectivity=26, mask=mask)
    
    plt.imshow(imageData[48], cmap='gray')
    plt.show()
    plt.imshow(labels[48])
    plt.show()
    
    return labels

def gaussian_smooth(myMolecule, sigma = 1):
    newMolecule = copy(myMolecule)
    smoothed = gaussian(newMolecule.data(), sigma=sigma)
    newMolecule.set_data(smoothed)
    return newMolecule
