import sys
import pathlib


sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
from general_utils.emd_utils import get_all_emd_name, get_associated_pdb

import matplotlib

from general_utils.graph_utils import remove_node_by_name, remove_edge_nodes
from general_utils.list_utils import combinations_i2jInK
from general_utils.workspace_utils import is_work_in_cluster
from process_graph.process_graph_utils import draw_graph_similarity, draw_graph_similarity_same_image

import matplotlib.pyplot as plt
from general_utils.database_utils import get_graph_pdb_db, get_chains_pdb_db, exists_mongo_db, get_all_archive_pdb, \
  memory_use, save_collection, load_collection, clear_collection
from process_graph.graph_algorithm import graph_aligning
import networkx as nx


import os

from general_utils.download_utils import download_pdb, download_emd_xml
from general_utils.pdb_utils import get_pdb_chain_sequence, get_similar_pdb_struct, get_similar_pdb_chain_structural, \
  get_similar_pdb_chain_sequential, get_chains_pdb, get_pdb_no_work
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains

# pdb = '5T4P'
# path = "./"

# pdb = '6VM4'

# download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
# chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))


# download_pdb(pdb, pdb_path)
# secuence = get_pdb_chain_sequence(pdb_path, "A")
# print(secuence)

# Todos se me parecen
# result_struct = get_similar_pdb_struct(pdb, -1)
# print(result_struct)


# result_chain_struct = get_similar_pdb_chain_structural(pdb, "A")
# print(result_chain_struct)
#
# for i in result_chain_struct:
#   pdb = i[0]
#
#   download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
#   chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))
#   pdb_to_mrc_chains(True, False, 5.0, '{0}/{1}.pdb'.format(path, pdb), path, chains,
#                     len(chains))
#   os.remove('{0}/{1}.pdb'.format(path, pdb))


# result_chain_sequence = get_similar_pdb_chain_sequential(pdb, "A")
# print(result_chain_sequence)

#
# import mrcfile
# import numpy as np
#
#
# def simulate_contour_level_value(mrc_file):
#   mrc = mrcfile.open(mrc_file)
#   data = mrc.data.flatten()
#
#   avg = np.average(data)
#   std = np.std(data)
#   min = np.min(data)
#   max = np.max(data)
#   var = np.var(data)
#   median = np.median(data)
#
#   rms = lambda V, axis=None: np.sqrt(np.mean(np.square(V), axis))
#   rms_value = rms(data)
#
#   multipli_value = (1 - (np.abs(min) / max)) * 2
#
#   try_value = (avg / 2) + \
#               median + \
#               std + \
#               var + \
#               multipli_value * rms_value
#
#   return try_value
#
# mrc = mrcfile.open('/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/main/6vof/6vof.mrc')
# #mrc = mrcfile.open('/home/lcastillo98/Desktop/emd_9334.map')
# #mrc = mrcfile.open('/home/lcastillo98/Desktop/emd_23065.map')
# #mrc = mrcfile.open('/home/lcastillo98/Desktop/emd_22776.map')
# data = mrc.data.flatten()
# #data = data[np.where( data > 0) ]
# print("Avg", np.average(data))
# print("STD", np.std(data))
# print("min", np.min(data))
# print("max", np.max(data))
# print("var", np.var(data))
# sumt = np.sum(data)
# print("sum", sumt)
# print("median", np.median(data))
# semi_sum = np.sum(data[np.where( data > 0.164 )])
# print("semisum", semi_sum)
#
# print("porcent masa", semi_sum/sumt)
#
# rms = lambda V, axis=None: np.sqrt(np.mean(np.square(V), axis))
# print("rms", rms(data))
#
# multipli_value = (1 - np.abs(np.min(data))/np.max(data))*2
# print(multipli_value)
#
# try_value = (np.average(data)/2) + \
#             np.median(data) + \
#             np.std(data) + \
#             np.var(data) + \
#             multipli_value*rms(data)
# print("try", try_value)
#
# 23065
# "Mas pequeno"
#
# print(simulate_contour_level_value('/home/lcastillo98/Desktop/EMD-5017.map'))
# print(get_pdb_no_work())

# chains = get_chains_pdb_db('1c5f')
# nums =[ 1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16]
# all = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
# test =[          'C', 'D', 'E',      'G', 'H', 'I',      'K', 'L', 'M',      'O'     ]
#
#
graph1 = get_graph_pdb_db('5w66', 10)
print(list(graph1.nodes))
#graph1 = get_graph_pdb_db('4c4g', 10)

graph2 = graph1.copy()

# remove_node_by_name(graph2, 13)
# remove_node_by_name(graph2, 15)
# remove_node_by_name(graph2, 17)
# remove_node_by_name(graph2, 3)
# remove_node_by_name(graph2, 4)
# remove_node_by_name(graph2, 6)



alignment_note, result = graph_aligning(graph1, graph2, 2, False)
print(alignment_note, result)

draw_graph_similarity_same_image(graph1, graph2, result)

#save_collection('./collection.json')
#clear_collection()
#load_collection('./collection.json')
#print(get_all_archive_pdb())
#
# result = combinations_i2jInK(10, 2, 5)
# print(result)
# print(len(result))

# print(is_work_in_cluster())
# print(exists_mongo_db())
# print(len(get_all_archive_pdb()))
# print(len(list(set(get_all_archive_pdb()))))
# print(memory_use())
# all_emd_names = get_all_emd_name()
# star_pos = all_emd_names.index('0561')
# for i in all_emd_names[star_pos:]:
#   result = get_associated_pdb(i)
#   if len(result)>1:
#     print(i)
#     print(result)
#     print("\n\n\n")
#   else:
#     print("no " + i)
#
#
# import sys
# import pathlib
#
#
# pathlib.Path(__file__).parent.absolute()
#
# sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../../em")
# import visualizer
# import processing
# import molecule
# import configparser
# import argparse
#
# import numpy as np
#
#
#
# myMolecule = molecule.Molecule("/tmp/tmpk4uj2qbj/1h1k/1h1k.mrc", recommendedContour=5)
# # Segment EM map with parameters steps (3) and sigma (1)
# myMolecule.generateSegments(3, 1)
# # Get generated dictionary, with contour ratio as keys and composited numpy array with segment masks and labels
# segments_at_contour_dict = myMolecule.getSegmentsMasks()
# # To get density array for each segment, we can use mask array indexing
# # What contour ratios do we have?
# # print(segments_at_contour_dict.keys())
# # Get segments masks at default contour ratio (which is 1)
# segments_masks = segments_at_contour_dict[1]
# # Print segment labels
# # print (segments_masks.keys())
# # Lets create a list to store a copy of densities for each segment
# result = []
# import utils._zernike as z
# zd = z.computeDescriptors(myMolecule.getEmMap().data())



# # Import math Library
# import math
#
# # Initialize the number of items to choose from
# n = 20
#
# # Print total number of possible combinations
#
# total = 0
#
# for i in range(n):
#   result = math.comb(n, i)
#   print(i, ":", result)
#   total += result
# print(total)