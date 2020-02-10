import numpy as np
from skimage.measure import regionprops
from skimage.transform import resize

class Model():

    defaultValue = 0
    indicatorValue = 1000


    def __init__(self,molecule, recommendedContour, cutoffRatios=[1]):
        self.molecule=molecule
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= len(cutoffRatios)

        moleculeData = self.molecule.data()
        dimentions = tuple([round(num) for num in self.getCellDim()])
        #c_map_array = np.ndarray((self.contoursNum,*dimentions))
        #map_data = resize(moleculeData,dimentions,order=3)
        c_map_array = np.ndarray((self.contoursNum,*self.getGridSize()))
        map_data = moleculeData
        for i,cutoffRatio in enumerate(self.cutoffRatios):
            map_at_contour = np.copy(map_data)
            map_at_contour[map_at_contour<cutoffRatio*self.contourLvl]=self.defaultValue
            map_at_contour[map_at_contour>=cutoffRatio*self.contourLvl]=self.indicatorValue
            c_map_array[i,:] = map_at_contour
        
        self.data_array = c_map_array
        self.segments=None
        self.atoms=None
        self.zDescriptors=None

    def getData(self):
        return self.data_array   

    def getCutoffLevels(self):
        return [cutoffRatio*self.contourLvl for cutoffRatio in self.cutoffRatios]

    def getGridSize(self):
        return self.molecule.grid_size()

    def getCellDim(self):
        return self.molecule.cell_dim()

    def getSegments(self):
        return self.labels,self.segments

    def getZernikeDescriptors(self, number):
        self.zDescriptors = np.ndarray(((self.contoursNum, number)))
        for i in range(self.contoursNum):
            if self.data[i].flags['C_CONTIGUOUS']:
                pass
            else:
                self.data[i,:] = self.data[i].ascontiguousarray('C') 
                



    #def setSegmentLabels(self, labels):
    #    regions = regionprops(labels,this.)