import numpy as np
from skimage.measure import regionprops
from skimage.transform import resize

from reader import Reader
import processing

## This module represents the molecule object with its properties and different data representations for each contour level. 

class Molecule():

    defaultValue = 0
    indicatorValue = 1000

    ## Initialize molecule object with a filename, recomended contour value, and a list of cut-off ratios.
    def __init__(self, filename, recommendedContour=None, cutoffRatios=[1]):
        ## Call reader module, a exception is expected to reaise if filename is not valid. 
        try:
            molecule_map = Reader.open(filename)
        except IOError:
            print("[ERROR]: Could not create Molecule object")
            exit
        map_data = molecule_map.data()
        contoursNum= len(cutoffRatios)
        dimentions = tuple([round(num) for num in map_data.getCellDim()])
        contour_maks = np.ndarray((contoursNum,*map_data.getGridSize()))

        for i,cutoffRatio in enumerate(self.cutoffRatios):
            data_at_contour = np.copy(map_data)
            data_at_contour[data_at_contour<cutoffRatio*self.contourLvl]=self.defaultValue
            data_at_contour[data_at_contour>=cutoffRatio*self.contourLvl]=self.indicatorValue
            contour_maks[i,:] = data_at_contour
        
        self.emMap=molecule_map
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= contoursNum

        self.contour_maks = contour_maks
        self.segments=None
        self.atoms=None
        self.zDescriptors=None


    def generateSegments(self, steps, sigma):
        #for range(self.contoursNum):
            
        #segments = processing.segment(steps, sigma)
        return 1

    def getContourMasks(self):
        return self.map_contours_data   

    def getCutoffLevels(self):
        return [cutoffRatio*self.contourLvl for cutoffRatio in self.cutoffRatios]

    def getGridSize(self):
        return self.map_data.grid_size()

    def getCellDim(self):
        return self.map_data.cell_dim()

    def getSegments(self):
        return self.segments


    def getZernikeDescriptors(self):
        self.zDescriptors = np.ndarray(((self.contoursNum, zernike_vec_len)))
        for i in range(self.contoursNum):
            if self.data[i].flags['C_CONTIGUOUS']:
                pass
            else:
                self.data[i,:] = self.data[i].ascontiguousarray('C') 
                
