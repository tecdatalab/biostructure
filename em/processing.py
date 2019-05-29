from scipy import ndimage as ndi
from skimage.morphology import watershed, erosion, dilation, ball
from skimage.feature import peak_local_max
from skimage.filters import  gaussian
from copy import copy
import numpy as np



def watershed_segmentation(myMolecule, level):
    imageData = myMolecule.data()
    thresholded = imageData > level
    kernel = np.ones((3,3,3))
    mask = erosion(thresholded)
    mask = dilation(thresholded)
    #-------
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    #-------
    distance = ndi.distance_transform_edt(thresholded)
    distance =  distance/np.linalg.norm(distance.ravel(), np.inf, keepdims=True)
    markers = ndi.label(np.where(distance>0.7, 1.0, 0.0), structure=kernel)[0]
    labels = watershed(-imageData, markers, mask=mask)
    '''
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xx, yy = np.mgrid[0:imageData.shape[2], 0:imageData.shape[1]]
    ax.plot_surface(xx, yy, -thresholded[48] ,rstride=1, cstride=1, cmap=plt.cm.gray, linewidth=0)
    plt.show()
    plt.imshow(imageData[48], cmap='gray')
    plt.show()
    plt.imshow(markers[48], cmap='gray')
    plt.show()
    plt.imshow(distance[48], cmap='gray')
    plt.show()
    plt.imshow(-distance[48], cmap='gray')
    plt.show()
    plt.imshow(labels[48])
    plt.show()
    '''
    return labels

def gaussian_smooth(myMolecule, sigma = 1):
    newMolecule = copy(myMolecule)
    smoothed = gaussian(newMolecule.data(), sigma=sigma)
    newMolecule.set_data(smoothed)
    return newMolecule
