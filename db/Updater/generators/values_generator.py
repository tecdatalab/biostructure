'''
Created on 28 abr. 2019
@author: luis98

Last modification on 7 aug. 2020
@author: dnnxl
'''
import sys
sys.path.append('../')

import os
import os.path
import requests
import numpy as np
from clint.textui import progress
from generators.molecule import Molecule
from generators.writer import Writer
from utilities import utility
from classes.segment_entry import Segment_entry, List_segments

dir = ""

def download_file(emd_id):
    URL = "http://ftp.wwpdb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz".format(emd_id)
    response = requests.get(URL, stream=True)
    with open('{0}temp/emd_{1}.map.gz'.format(dir, emd_id), 'wb') as file:
        total_length = int(response.headers.get('content-length'))
        for chunk in progress.bar(response.iter_content(chunk_size = 8192), expected_size=(total_length/8192) + 1):
            if chunk:
                file.write(chunk)
                file.flush()
                os.fsync(file.fileno())
    os.system("gzip -d -f {0}temp/emd_{1}.map.gz".format(dir, emd_id))

def get_min_max_density(emd_id, contour):
    text = os.popen(
        './{0}em_volume {1} {2}temp/emd_{3}.map'.format(dir, contour, dir, emd_id)).read()
    list_a = text.split()
    return (float(list_a[10]), float(list_a[11]))

def generate_descriptor_file(emd_id, contour, exit_file_name):
    execute_string = "./{0}map2zernike {0}temp/emd_{1}.map -n 20 -c {2} -p {0}temp/{3} ".format(
            dir,
            emd_id,
            contour,
            exit_file_name)
    os.system(
        execute_string)
    if (not os.path.exists("{0}temp/{1}.inv".format(dir, exit_file_name))):
        raise Exception('Emd descriptor file was not created.')

def generate_descriptor_file_for_union(emd_id, contour_list, exit_file_name):
    temp_file_name = "temp"
    for i in contour_list:
        generate_descriptor_file(emd_id, i, temp_file_name)
        os.system(
            "tail -$(wc -l {0}temp/{1}.inv | awk '{{print $1-1}}') {0}temp/{1}.inv >> {0}temp/{2}.inv ".format(
                dir,
                temp_file_name,
                exit_file_name))

    os.system("rm {0}temp/{1}.inv ".format(dir, temp_file_name))

def generate_descriptors_files(emd_id, contour, std):
    max_density = get_min_max_density(emd_id, contour)[1]
    one_third_contour = contour + ((max_density - contour) / 3)
    two_thirds_contour = contour + ((max_density - contour) * (2 / 3))
    one_std_contour = contour + std

    generate_descriptor_file(emd_id, contour, "contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, one_third_contour], "one_third_contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, two_thirds_contour], "two_thirds_contour")
    generate_descriptor_file_for_union(
        emd_id, [
            contour, one_third_contour, two_thirds_contour], "one_third_two_thirds_contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, one_std_contour], "one_std_contour")

def remove_map(emd_id):
    os.system("rm -rf {0}temp ".format(dir))
    os.system("mkdir {0}temp ".format(dir))

def get_emd_descriptors(emd_id, contour, std):

    files = [
        "contour",
        "one_third_contour",
        "two_thirds_contour",
        "one_third_two_thirds_contour",
        "one_std_contour"]
    result = []

    generate_descriptors_files(emd_id, contour, std)

    for i in files:
        with open("{0}temp/{1}.inv".format(dir, i)) as f:
            array = [float(line) for line in f]
            result.append(array)

    for i in files:
        os.system("rm {0}temp/{1}.inv ".format(dir, i))

    return result

def get_volume_map(emd_id, recommendedContour, cutoffRatios):
    '''
        Returns the list of strings of the volume map, following
            the next format example : 
                ["{contour : volume}", ..., "{contour_n : volume_n}"]
        Parameters:
            emd_id (int)                 -  The emd id.
            recommendedContour (int)     -  The recommended contour level of the map.
            cutoffRatios(list of float)  -  The cut off ratios of the map.
        Returns:
            map_volume (list of strings)  -  The list of strings of contour with the 
                                             respective volume values.
    '''
    temp_molecule = Molecule("{0}temp/emd_{1}.map".format(dir, emd_id), recommendedContour=recommendedContour, cutoffRatios=cutoffRatios)
    map_volume = utility.format_volume_map(temp_molecule.getVolume())
    return map_volume

def get_molecule(emd_id, recommendedContour, cutoffRatios):
    '''
        Returns the a Molecule object with the emd_id and 
            recommended contour and cut off ratios.
        Parameters:
            emd_id (int)                 -  The emd id.
            recommendedContour (int)     -  The recommended contour level of the map.
            cutoffRatios(list of float)  -  The cut off ratios of the map.
        Returns:
            molecule (object)  -  The molecule object.
    '''
    molecule = Molecule('{0}temp/emd_{1}.map'.format(dir, emd_id), recommendedContour=recommendedContour, cutoffRatios=cutoffRatios)
    return molecule

def generate_segments(molecule, emd_id, map_id, algorithm_id, contour=1):
    '''
        Returns an object List_segments with the segments generated 
            by the default algortihm. 
        Parameters:
            molecule (object)  -  The object molecule.
            emd_id (int)       -  The emd id.
            map_id (int)       -  The map id.
            algorithm_id(int)  -  The algorithm id.
            contour(int)       -  The contour.
        Returns:
            List_segments (object)  -  The List_segments object.
    '''
    # Segment EM map with parameters steps (3) and sigma (1).
    molecule.generateSegments(3,1)
    # Get generated dictionary, with contour ratio as keys and composited numpy array with segment masks and labels.
    segments_at_contour_dict = molecule.getSegmentsMasks()
    # Get segments masks at default contour ratio (which is 1).    
    segments_masks = segments_at_contour_dict[contour]
    # Lets create a list to store a copy of densities for each segment.    
    segment_dict = {}
    # Get the volume of the segments with the respective contour.
    segments_volume_at_contour_dict = molecule.getSegmentsVolume()
    # Get segments volumes at default contour ratio (which is 1)    
    segments_volume = segments_volume_at_contour_dict[contour]
    # Initialize a writer object to write the maps.
    writer = Writer()
    segments = []
    for key_mask, key_volume in zip(segments_masks, segments_volume):
        # Create a copy of map densities
        densities = np.copy(molecule.getEmMap().data())
        # Set voxels outside segment to 0 
        densities[np.logical_not(segments_masks[key_mask])] = 0
        segment_dict[key_mask] = densities
        # The directory where the segments generate will be stored.
        directory = '{0}segments/{1}/{2}'.format(dir, emd_id, contour)
        # Created the directory with the format segments/emd_id/contour/
        os.makedirs(directory, exist_ok=True)
        # The segment file name with the convention seg_id.map
        segment_file_name = '{0}/seg_{1}.map'.format(directory, key_mask)
        # Write the segment in the file.
        writer.write(segment_file_name, molecule.emMap, densities.T)
        # Get the path segment.
        path_segment = os.path.dirname(os.path.abspath(segment_file_name))
        # Create the object segment and append to the list of segments.
        seg_temp = Segment_entry(
            None, map_id, algorithm_id, segments_volume[key_volume], contour, path_segment)
        segments.append(seg_temp)
        
    return List_segments(segments)


#download_file(1010)
#emd_id = 1010
#molecule = get_molecule(emd_id, 7, [1,0.5])
#generate_segments(molecule, emd_id)
"""
def segmentation_default(emd_id, recommendedContour, cutoffRatios):
    temp_molecule = Molecule("{0}temp/emd_{1}.map".format(dir, emd_id), recommendedContour=recommendedContour, cutoffRatios=cutoffRatios)
    temp_molecule.generateSegments(3,1)
    temp_segment_volumes = temp_molecule.getSegmentsVolume()
    return temp_segment_volumes
"""
# download_file("0001")
# print(get_min_max_density("0001",0.018))
# generate_descriptor_file("0001",0.018,"normal_contour")
#result = get_emd_descriptors("0001",0.018, 0.0077065425)

# for i in result:
#    print(i)