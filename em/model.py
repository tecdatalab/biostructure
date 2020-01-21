import numpy as np

class Model():
    def __init__(self,molecule, recommendedContour, cutoffRatios=[1], segments=None, atoms=None):
        self.molecule=molecule
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= len(cutoffRatios)
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

