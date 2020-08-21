import sys
import pathlib

pathlib.Path(__file__).parent.absolute()
from process_mrc.segment import Segment
import glob

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../../em")
import visualizer
import processing
import molecule
import configparser
import argparse

import numpy as np


def get_mrc_segments(mrc_path, recommendedContour_p, steps, sigma):
    # Initialize molecule object with arguments: filename, recomended contour value and an optional list of cut-off
    # ratios.
    myMolecule = molecule.Molecule(mrc_path, recommendedContour=recommendedContour_p)
    # Segment EM map with parameters steps (3) and sigma (1)
    myMolecule.generateSegments(steps, sigma)
    # Get generated dictionary, with contour ratio as keys and composited numpy array with segment masks and labels
    segments_at_contour_dict = myMolecule.getSegmentsMasks()
    # To get density array for each segment, we can use mask array indexing
    # What contour ratios do we have?
    # print(segments_at_contour_dict.keys())
    # Get segments masks at default contour ratio (which is 1)
    segments_masks = segments_at_contour_dict[1]
    # Print segment labels
    # print (segments_masks.keys())
    # Lets create a list to store a copy of densities for each segment
    result = []
    for key in segments_masks:
        # Create a copy of map densities
        densities = np.copy(myMolecule.getEmMap().data())
        # Set voxels outside segment to 0
        densities[np.logical_not(segments_masks[key])] = 0
        result.append(Segment(key, densities, None))

    # then you can compute zernike descriptors for each segment, lets create a dict to store descriptors for each
    # segment # lets import the module
    import utils._zernike as z
    # lets create a dictionary to store descriptors for each segment
    for i in result:
        zd = z.computeDescriptors(i.mask)
        i.zd_descriptors = zd

    # print("Can_points", len(np.where(myMolecule.getEmMap().data() > 0)[0]))

    return result, myMolecule.getEmMap().data().shape


def get_mrc_synthetic_segments_pdb(folder_segments, recommendedContour_p):
    segments_paths = list(glob.glob(folder_segments + "/*_*.mrc"))
    # print(segments_paths)
    actual_id = 1
    result = []
    len_X = 0
    len_Y = 0
    len_Z = 0
    for path in segments_paths:
        # print(path)
        ## Initialize molecule object with arguments: filename, recomended contour value and an optional list of cut-off ratios.
        myMolecule = molecule.Molecule(path, recommendedContour=recommendedContour_p)
        densitie = np.copy(myMolecule.getEmMap().data())
        result.append(Segment(actual_id, densitie, None))
        actual_id += 1
        len_X = max(len_X, myMolecule.getEmMap().data().shape[0])
        len_Y = max(len_Y, myMolecule.getEmMap().data().shape[1])
        len_Z = max(len_Z, myMolecule.getEmMap().data().shape[2])
        # print(myMolecule.contour_masks.shape)

    for i in result:
        zd = z.computeDescriptors(i.mask)
        i.zd_descriptors = zd

    return result, (len_X, len_Y, len_Z)


def get_mrc_one(mrc_path, recommendedContour_p):
    # Initialize molecule object with arguments: filename, recomended contour value and an optional list of cut-off
    # ratios.
    actual_id = 1
    result = []
    myMolecule = molecule.Molecule(mrc_path, recommendedContour=recommendedContour_p)
    densities = np.copy(myMolecule.getEmMap().data())
    # Set voxels outside segment to 0
    densitie = np.copy(myMolecule.getEmMap().data())
    result.append(Segment(actual_id, densitie, None))

    # then you can compute zernike descriptors for each segment, lets create a dict to store descriptors for each
    # segment # lets import the module
    import utils._zernike as z
    # lets create a dictionary to store descriptors for each segment
    for i in result:
        zd = z.computeDescriptors(i.mask)
        i.zd_descriptors = zd

    # print("Can_points", len(np.where(myMolecule.getEmMap().data() > 0)[0]))

    return result, myMolecule.getEmMap().data().shape

# result = get_mrc_segments("../pdb_mrc/exit_pdb/175d/175d.mrc", 7, 3, 1)
# result = get_mrc_synthetic_segments_pdb("../pdb_mrc/exit_pdb/175d", 7)
