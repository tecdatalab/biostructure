import numpy as np
from skimage.measure import regionprops
from skimage.transform import resize

class Model():
    def __init__(self,molecule, recommendedContour, cutoffRatios=[1]):
        self.molecule=molecule
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= len(cutoffRatios)
        self.segments=None
        self.data=None
        self.atoms=None
        self.zDescriptors=None

    def getData(self):
        moleculeData = self.molecule.data()
        dimentions = tuple([round(num) for num in self.getCellDim()])
        mydata = np.ndarray((self.contoursNum,*dimentions))
        
        for i,cutoffRatio in enumerate(self.cutoffRatios):
            data = resize(moleculeData,dimentions,order=3)
            data[data<cutoffRatio*self.contourLvl]=0
            mydata[i,:] = data
        self.data = mydata
        return self.data   

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