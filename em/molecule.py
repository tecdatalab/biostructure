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
            reader = Reader()
            reader.open(filename)
            molecule_map = reader.read()
        except IOError:
            print("[ERROR]: Could not create Molecule object")
            exit()
        map_data = molecule_map.data()
        contoursNum= len(cutoffRatios)
        contour_masks = np.ndarray((contoursNum,*molecule_map.grid_size()))
        
        self.emMap=molecule_map
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= contoursNum

        for i,cutoffRatio in enumerate(self.cutoffRatios):
            data_at_contour = np.copy(map_data)
            data_at_contour[data_at_contour<cutoffRatio*self.contourLvl]=self.defaultValue
            data_at_contour[data_at_contour>=cutoffRatio*self.contourLvl]=self.indicatorValue
            contour_masks[i,:] = data_at_contour

        self.contour_masks = contour_masks
        self.seg_masks=None
        self.labels=None
        self.atoms=None
        self.zDescriptors=None


    def generateSegments(self, steps, sigma):
        labels = processing.segment(self.emMap, steps, sigma, self.getCutoffLevels(), self.contoursNum)
        volume_reg_lvl = []
        voxel_vol = self.emMap.voxelVol()
        levels_masks = []
        for lvl_labels in labels:
            volume_reg_dict = {}
            regprops = regionprops(lvl_labels)
            segment_masks = np.zeros(len(regprops), [("label", np.int, 1), ("mask", np.int, self.emMap.grid_size())])
            for i,reg in enumerate(regprops):
                volume_reg_dict[reg.label] = voxel_vol*reg.area
                segment_masks["label"][i]=reg.label
                segment_indexes=reg.coords
                segment_masks["mask"][i][segment_indexes[:,0],segment_indexes[:,1],segment_indexes[:,2]] = 1
            volume_reg_lvl.append(volume_reg_dict)
            levels_masks.append(segment_masks)
        self.labels = labels
        self.seg_masks = levels_masks
        print(self.emMap.stdev())



    def getContourMasks(self):
        return self.contour_masks   

    def getEmMap(self):
        return self.emMap

    def getCutoffLevels(self):
        return [cutoffRatio*self.contourLvl for cutoffRatio in self.cutoffRatios]

    def getGridSize(self):
        return self.emMap.grid_size()

    def getCellDim(self):
        return self.emMap.cell_dim()

    def getSegments(self):
        return self.segments


    def getZernikeDescriptors(self):
        self.zDescriptors = np.ndarray(((self.contoursNum, zernike_vec_len)))
        for i in range(self.contoursNum):
            if self.contour_maks[i].flags['C_CONTIGUOUS']:
                pass
            else:
                self.contour_maks[i,:] = self.data[i].ascontiguousarray('C') 
                
