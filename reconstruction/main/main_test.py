
import re
import sys
import pathlib


sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
from pymol import cmd
from general_utils.temp_utils import clean_work_dir

from experiment.utils_general import pdb_percentage, check_RMSD_result_algorithm, make_dir_pdb
from general_utils.cif_utils import get_chains_cif, get_cif_chain_sequence
from to_mrc.cif_2_mrc import cif_to_mrc_chains

from general_utils.emd_utils import get_all_emd_name, get_associated_pdb

import matplotlib

from general_utils.graph_utils import remove_node_by_name, remove_edge_nodes
from general_utils.list_utils import combinations_i2jInK, seudo_combinations_i2jInK
from general_utils.workspace_utils import is_work_in_cluster
from process_graph.process_graph_utils import draw_graph_similarity, draw_graph_similarity_same_image

import matplotlib.pyplot as plt
from general_utils.database_utils import get_graph_pdb_db, get_chains_pdb_db, exists_mongo_db, get_all_archive_pdb, \
  memory_use, save_collection, load_collection, clear_collection, get_chain_to_number_chain, delete_pdb_db, \
  get_online_sequences, get_sequence_pdb_db, get_online_chain_map_pdb2cif, get_pdb2cif_db
from process_graph.graph_algorithm import graph_aligning
import networkx as nx

import os

from general_utils.download_utils import download_pdb, download_emd_xml
from general_utils.pdb_utils import get_pdb_chain_sequence, get_similar_pdb_struct, get_similar_pdb_chain_structural, \
  get_similar_pdb_chain_sequential, get_chains_pdb, get_pdb_no_work
from to_mrc.pdb_2_mrc import pdb_to_mrc_chains
import general_utils

if is_work_in_cluster():
  general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_main_test"
else:
  general_utils.temp_utils.global_temp_dir = "./main_test"
clean_work_dir()

# pdb = '5T4P'

#
# download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
# chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))
#
#
# download_pdb(pdb, pdb_path)
# secuence = get_pdb_chain_sequence(pdb_path, "A")
# print(secuence)
#
# #Todos se me parecen
# result_struct = get_similar_pdb_struct(pdb, -1)
# print(result_struct)
#
#
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
# graph1 = get_graph_pdb_db('5w5y', 10)
# print(list(graph1.nodes))
#
# graph2 = graph1.copy()
#
# # remove_node_by_name(graph2, 13)
# # remove_node_by_name(graph2, 15)
# # remove_node_by_name(graph2, 17)
# # remove_node_by_name(graph2, 3)
# # remove_node_by_name(graph2, 4)
# # remove_node_by_name(graph2, 6)
#
#
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
#
# draw_graph_similarity_same_image(graph1, graph2, result)

# save_collection('./collection.json')
# clear_collection()
# load_collection('./collection.json')
# print(get_all_archive_pdb())
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
# total = combinations_i2jInK(1000, 1, 5, check_max=True)
#
# print(total)


# Gen PDB with chains

# path = "./"
#
# pdb = '7no3'
#
# download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
# chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))
#
# pdb_to_mrc_chains(True, False, 5.0, '{0}/{1}.pdb'.format(path, pdb), path, chains,
#                     len(chains))


'''
1is8
2HOD
3LEL
6h7w
'''

'''
1is8

Remove
P, Q, R, S, T

'''

# pdb = '1is8'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'P'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'Q'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'R'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'S'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'T'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)


'''
1is8

Remove
K, L, M, N, O

'''

# pdb = '1is8'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'K'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'L'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'M'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'N'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'O'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)


'''
1is8

Remove
P, Q, R, S, T, K, L, M, N, O

'''
#
# pdb = '1is8'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'P'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'Q'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'R'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'S'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'T'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'K'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'L'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'M'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'N'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'O'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)


'''
2HOD
'''

'''
2hod

Remove
A, B, C, D, E, F, M, N, O, P

'''

# pdb = '2hod'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'A'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'B'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'C'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'D'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'E'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'F'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'M'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'N'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'O'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'P'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)

'''
2hod

Remove
G, H, I, J, K, L, Q, R, S, T

'''

# pdb = '2hod'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'G'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'H'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'I'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'J'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'K'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'L'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'Q'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'R'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'S'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'T'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)


'''
3lel
'''

'''
3lel

Remove
A, B, C, D, E, F, G, H, I, J

'''

# pdb = '3lel'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'A'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'B'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'C'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'D'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'E'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'F'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'G'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'H'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'I'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'J'))
#
# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)


'''
3lel

Remove
K, L, M, N, O, P, Q, R, S, T

'''

# pdb = '3lel'
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'A'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'B'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'C'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'D'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'E'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'F'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'G'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'H'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'I'))
# remove_node_by_name(graph2, get_chain_to_number_chain(pdb, 'J'))
#
# load_collection("/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/main/collection.json")

# alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note, result)
# draw_graph_similarity_same_image(graph1, graph2, result)
# pdb = "5taw"
# # download_pdb(pdb, "./7ndg.pdb")
# # result = get_chains_pdb("./7ndg.pdb")
# # print(result)
# # delete_pdb_db(pdb)
# graph1 = get_graph_pdb_db(pdb, 4)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
# for i in get_chains_pdb_db(pdb):
#   print("Chain=", i)
#   result = get_similar_pdb_chain_sequential(pdb, i)
#   print(result)

# path = "./"
# #
# pdb = '4v4r'
# #
# download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
# chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))
# print(chains)
# #
# # pdb_to_mrc_chains(True, False, 5.0, '{0}/{1}.pdb'.format(path, pdb), path, chains,
# #                     len(chains))
#
# pdb = '4v4r'
# graph1 = get_graph_pdb_db(pdb, 8)
# print(list(graph1.nodes))
# print(get_chains_pdb_db(pdb))
#
# graph2 = graph1.copy()
#
# # alignment_note, result = graph_aligning(graph1, graph2, 2, False)
# # print(alignment_note, result)
# # draw_graph_similarity_same_image(graph1, graph2, result)
#
# print(graph1.nodes[1])
# print(len(get_all_archive_pdb()))
# cif_path = "./4v4r.cif"
# cif_path = "./3j79.cif"
# pdb_path = "./3j79.cif"
# chains = get_chains_cif(cif_path)
# print(chains)
# #
# cif_to_mrc_chains(True, False, 5.0, cif_path, "./", chains, len(chains))
#
# sequence = get_cif_chain_sequence(cif_path, "4v4r", "A")
# print(sequence)
# print("\n\n\n")

# sequence = get_sequence_pdb_db("1brs", "A")
# print(sequence)

# all_pdb_names = pdb_percentage(10)

# pdb_name = "3whe"
# chains = get_chains_pdb_db(pdb_name)
# print(chains)
# result = get_online_sequences(pdb_name)
# #
# for i in result.keys():
#   print(i)
#   print(result[i])
# chains_sec = "AB[B,C],C[D],E[F],G,P"
# chains_sec = "BA[B,C],C[D],E[F],G,P"
#
# #chains = re.split(r'\],|,', chains_sec)
# chains = re.findall(r'[A-Za-z0-9]+\[[^]]*\]|[A-Za-z0-9]+', chains_sec)
#
# chain_cif = re.findall(r'\[[^]]*\]', chains[0])
#
# print(chains)
# print(chain_cif)


# pdb_name = "1vy7"
# # delete_pdb_db(pdb_name)
# # chains = get_chains_pdb_db(pdb_name)
# # print(chains)
# equivalent = get_pdb2cif_db(pdb_name)
# print(equivalent)
# # result = get_online_sequences(pdb_name)
#
# result = get_pdb2cif_db("7lj8")
# print(result)

# all_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R']
#
# check_list = [
#   [2, 15],
#   [1, 1],
#   [6, 6],
#   [16, 16],
#   [13, 13],
#   [4, 3],
#   [12, 5],
#   [11, 10],
#   [5, 9],
#   [8, 2],
#   [14, 8],
#   [17, 14],
#   [10, 18],
#   [15, 11],
#   [18, 4],
#   [9, 12],
#   [7, 17]]
#
# original_pdb = "4u0d"
# changed_pdb = "1gwp"
# changed_chain = "G"
# pdb_work_chain = "A"
#
work_dir = "./RMSD"
# # #
# # # # result = check_RMSD_result_algorithm(work_dir, all_chains, check_list, original_pdb, changed_pdb,
# # # #                                      changed_chain, pdb_work_chain)
# # # # print(result)
# # #
make_dir_pdb(work_dir, "3sdd")
#cif_path = "./4ayb.cif"
#chains = get_chains_cif(cif_path)
#print(chains)




pdbf1 = '/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/main/RMSD/3sdd/3sdd_D.pdb';
pdbf2 = '/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/main/RMSD/3sdd/3sdd_D.pdb'

cmd.load(pdbf1, "3sddi_D")
cmd.load(pdbf2, "3sddj_D")

result = cmd.align("3sddi_D", "3sddj_D")
print(result)
