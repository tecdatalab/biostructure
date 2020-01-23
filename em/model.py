import numpy as np
from skimage import regionprops

class Model():
    def __init__(self,molecule, recommendedContour, cutoffRatios=[1], segments=None, atoms=None):
        self.molecule=molecule
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= len(cutoffRatios)
        self.data_
        self.segments=[]
        self.atoms=atoms

    def getData(self):
        moleculeData = self.molecule.data()
        mydata = np.ndarray((self.contoursNum,*moleculeData.shape))
        
        for i,cutoffRatio in enumerate(self.cutoffRatios):
            data = np.copy(moleculeData)
            data[data<cutoffRatio*self.contourLvl]=0
            mydata[i,:] = data
        return mydata   

    def getCutoffLevels(self):
        return [cutoffRatio*self.contourLvl for cutoffRatio in self.cutoffRatios]

    def getGridSize(self):
        return self.molecule.grid_size()

    def getCellDim(self):
        return self.molecule.cell_dim()

    def getSegments(self):
        return self.labels,self.segments

    #def setSegmentLabels(self, labels):
    #    regions = regionprops(labels,this.)