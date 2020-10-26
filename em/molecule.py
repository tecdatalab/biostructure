import numpy as np
from skimage.measure import regionprops
from skimage.transform import resize

from em.reader import Reader
from em import processing

## This module represents the molecule object with its properties and different data representations for each contour level. 

class Molecule():

    defaultValue = 0
    indicatorValue = 1

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
        # Use dictionary instead
        #contour_masks = np.ndarray((contoursNum,*molecule_map.grid_size()))
        contour_masks = {}
        self.emMap=molecule_map
        self.contourLvl=recommendedContour
        self.cutoffRatios=cutoffRatios
        self.contoursNum= contoursNum

        for cutoffRatio in self.cutoffRatios:
            lvl = cutoffRatio*self.contourLvl
            mask_at_contour = np.zeros(map_data.shape, dtype=np.bool)
            mask_at_contour[map_data>=lvl]=self.indicatorValue
            contour_masks[cutoffRatio] = mask_at_contour

        self.contour_masks = contour_masks
        self.seg_masks=None
        self.labels=None
        self.atoms=None
        self.zDescriptors=None

    # step and sigma are input arguments, set seg_masks which stores a collection of segment maks as composited numpy array for each countour level. 
    # Returns a dictionary where key is level ratio and value is the composited numpy array
    def generateSegments(self, steps, sigma):
        cutoffLvls = self.getCutoffLevels()
        labels = processing.segment(self.emMap, steps, sigma, cutoffLvls, self.contoursNum)
        volume_reg_lvl = []
        voxel_vol = self.emMap.voxelVol()
        molecule_masks = {}
        for lvl,labels in zip(self.cutoffRatios,labels):
            volume_reg_dict = {}
            level_masks = {}
            regprops = regionprops(labels)
            for i,reg in enumerate(regprops):
                segment_mask = np.zeros(self.emMap.grid_size(), dtype=np.int)
                segment_indexes=reg.coords
                segment_mask[segment_indexes[:,0],segment_indexes[:,1],segment_indexes[:,2]] = 1
                volume_reg_dict[reg.label] = voxel_vol*reg.area
                level_masks[reg.label]=segment_mask
            volume_reg_lvl.append(volume_reg_dict)
            molecule_masks[lvl]= level_masks
        self.labels = labels
        self.seg_masks = molecule_masks



    def getContourMasks(self):
        return self.contour_masks   

    def getEmMap(self):
        return self.emMap

    def getCutoffLevels(self):
        return {cutoffRatio: cutoffRatio*self.contourLvl for cutoffRatio in self.cutoffRatios}

    def getGridSize(self):
        return self.emMap.grid_size()

    def getCellDim(self):
        return self.emMap.cell_dim()

    def getSegmentsMasks(self):
        return self.seg_masks


    def getZernikeDescriptors(self):
        self.zDescriptors = np.ndarray(((self.contoursNum, zernike_vec_len)))
        for i in range(self.contoursNum):
            if self.contour_maks[i].flags['C_CONTIGUOUS']:
                pass
            else:
                self.contour_maks[i,:] = self.data[i].ascontiguousarray('C') 

    def getVolume(self):
        volume_contour_dict = dict()
        voxel_vol = self.emMap.voxelVol()
        for cutoffRatio in self.cutoffRatios:
            mask_at_level = self.contour_masks[cutoffRatio]
            volume_contour_dict[cutoffRatio] = np.sum(mask_at_level)*voxel_vol
        return volume_contour_dict

    def getSegmentsVolume(self):
        volume_contour_dict = dict()
        voxel_vol = self.emMap.voxelVol()
        segment_masks = self.getSegmentsMasks()
        for cutoffRatio in self.cutoffRatios:
            mask_at_level = segment_masks[cutoffRatio]
            volume_contour_dict[cutoffRatio] = {label_id : np.sum(mask_at_level[label_id])*voxel_vol for label_id in mask_at_level.keys()}
        return volume_contour_dict


    def getVoxelVol(self):
        return self.emMap.voxelVol()

    def getVoxelSize(self):
        return self.emMap.voxelSize()

    def getCorrelation(self, molecule, levels=[1]):
        if(self.getGridSize()!=molecule.getGridSize()):
            raise ValueError("Both maps must have same dimensions")
        else:
            result = {}
            moleculeA_masks = self.getContourMasks()
            moleculeB_masks = molecule.getContourMasks()
         
            for level in levels:
                mask_A = moleculeA_masks[level]
                mask_B = moleculeB_masks[level]
                data_A = np.multiply(self.emMap.data(), mask_A)
                data_B = np.multiply(self.emMap.data(), mask_B)
                mean_A = np.mean(data_A)
                mean_B = np.mean(data_B)
                corr_numerator = np.sum(np.einsum("ijk, ijk -> ijk", (data_A - mean_A), (data_B - mean_B)))
                corr_denominator =  np.sqrt( np.sum(np.power((data_A - mean_A), 2)) * np.sum(np.power((data_B - mean_B), 2)))
                result[level] = corr_numerator/corr_denominator
            return result
            
    def getOverlap(self, molecule, levels=[1]):
        if(self.getGridSize()!=molecule.getGridSize()):
            raise ValueError("Both maps must have same dimensions")
        else:
            result = {} 
            moleculeA_masks = self.getContourMasks()
            moleculeB_masks = molecule.getContourMasks()

            for level in levels:
                mask_A = moleculeA_masks[level]
                mask_B = moleculeB_masks[level]
                overlap = np.sum(np.einsum("ijk, ijk -> ijk", mask_A, mask_B))
                result[level] = overlap
            return result 
         



                
