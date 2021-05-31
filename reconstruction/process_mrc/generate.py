import sys
import pathlib

from general_utils.mrc_uilts import get_mrc_level

pathlib.Path(__file__).parent.absolute()
from process_mrc.segment import Segment
from process_mrc.biomolecular_structure import Biomolecular_structure
import glob

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../../em")
import visualizer
import processing
import molecule
import configparser
import argparse

import numpy as np


def get_mrc_segments(mrc_path, steps, sigma, recommendedContour_p=None, calculate_Z3D=True):
  if recommendedContour_p == None:
    recommendedContour_p = get_mrc_level(mrc_path)
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
  if calculate_Z3D:
    for i in result:
      zd = z.computeDescriptors(i.mask)
      i.zd_descriptors = zd

  # print("Can_points", len(np.where(myMolecule.getEmMap().data() > 0)[0]))
  original_structure = \
    Biomolecular_structure(myMolecule.getEmMap().data(), z.computeDescriptors(myMolecule.getEmMap().data()))

  return result, original_structure


def get_mrc_synthetic_segments_pdb(mrc_path, folder_segments, recommendedContour_p=None, calculate_Z3D=True):
  if recommendedContour_p == None:
    recommendedContour_p = get_mrc_level(mrc_path)
  myMolecule_complete = molecule.Molecule(mrc_path, recommendedContour=recommendedContour_p)
  segments_paths = list(glob.glob(folder_segments + "/*_*.mrc"))
  # print(segments_paths)
  actual_id = 1
  result = []
  for path in segments_paths:
    # print(path)
    ## Initialize molecule object with arguments: filename, recomended contour value and an optional list of cut-off ratios.
    myMolecule = molecule.Molecule(path, recommendedContour=recommendedContour_p)
    densitie = np.copy(myMolecule.getEmMap().data())
    result.append(Segment(actual_id, densitie, None))
    actual_id += 1
    # print(myMolecule.contour_masks.shape)

  import utils._zernike as z
  if calculate_Z3D:
    for i in result:
      zd = z.computeDescriptors(i.mask)
      i.zd_descriptors = zd

  original_structure = \
    Biomolecular_structure(myMolecule_complete.getEmMap().data(),
                           z.computeDescriptors(myMolecule_complete.getEmMap().data()))

  return result, original_structure


def get_mrc_synthetic_segments_pdb_list(mrc_path, folder_segments, list_segments, recommendedContour_p=None, calculate_Z3D=True):
  if recommendedContour_p == None:
    recommendedContour_p = get_mrc_level(mrc_path)
  myMolecule_complete = molecule.Molecule(mrc_path, recommendedContour=recommendedContour_p)
  segments_paths = [folder_segments + i for i in list_segments]
  # print(segments_paths)
  actual_id = 1
  result = []
  for path in segments_paths:
    # print(path)
    ## Initialize molecule object with arguments: filename, recomended contour value and an optional list of cut-off ratios.
    myMolecule = molecule.Molecule(path, recommendedContour=recommendedContour_p)
    densitie = np.copy(myMolecule.getEmMap().data())
    result.append(Segment(actual_id, densitie, None))
    actual_id += 1
    # print(myMolecule.contour_masks.shape)

  import utils._zernike as z
  if calculate_Z3D:
    for i in result:
      zd = z.computeDescriptors(i.mask)
      i.zd_descriptors = zd

    original_structure = \
      Biomolecular_structure(myMolecule_complete.getEmMap().data(),
                             z.computeDescriptors(myMolecule_complete.getEmMap().data()))

  return result, original_structure


def get_mrc_one(mrc_path, recommendedContour_p=None, calculate_Z3D=True):
  if recommendedContour_p == None:
    recommendedContour_p = get_mrc_level(mrc_path)
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
  if calculate_Z3D:
    for i in result:
      zd = z.computeDescriptors(i.mask)
      i.zd_descriptors = zd

  # print("Can_points", len(np.where(myMolecule.getEmMap().data() > 0)[0]))
  original_structure = \
    Biomolecular_structure(myMolecule.getEmMap().data(),
                           result[0].zd_descriptors)

  return result, original_structure
