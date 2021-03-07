import sys
import pathlib

import matplotlib

from process_graph.process_graph_utils import draw_graph_similarity

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import matplotlib.pyplot as plt
from general_utils.database_utils import get_graph_pdb_db, get_chains_pdb_db
from process_graph.graph_algorithm import graph_aligning
from process_graph.graph_algorithm2 import graph_aligning2
from process_graph.graph_algorithm3 import graph_aligning3
import networkx as nx


import os

from general_utils.download_utils import download_pdb
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

chains = get_chains_pdb_db('1c5f')
nums =[ 1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16]
all = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
test =[          'C', 'D', 'E',      'G', 'H', 'I',      'K', 'L', 'M',      'O'     ]


print(chains)
graph1 = get_graph_pdb_db('1c5f', 10)
graph2 = get_graph_pdb_db('1c5f', 10)



##Remove for test
graph2.remove_node(12) #ok
graph2.remove_node(15) #ok
graph2.remove_node(7)
graph2.remove_node(3)
#graph2.remove_node(1)
#graph2.remove_node(4)
#graph2.remove_node(5)
#graph2.remove_node(6)
#graph2.remove_node(7)
#graph2.remove_node(11)
#graph2.remove_node(12)
graph2.remove_edge(13, 6)
graph2.remove_edge(13, 5)
graph2.remove_edge(13, 11)
graph2.remove_edge(14, 11)


matplotlib.use('TKAgg')


# alignment_note1, result1 = graph_aligning(graph1, graph2, 2, False)
# print(alignment_note1, result1)
# print("\n\n\n\n")
# alignment_note2, result2 = graph_aligning2(graph1, graph2, 2, False)
# print(alignment_note2, result2)
# print("\n\n\n\n")
print(graph2.nodes())
alignment_note3, result3 = graph_aligning3(graph1, graph2, 2, False)
print(alignment_note3, result3)


plt.figure(1)
plt.suptitle('Graph 1', fontsize=16)
nx.draw(graph1, with_labels=True)
plt.figure(2)
plt.suptitle('Graph 2', fontsize=16)
nx.draw(graph2, with_labels=True)
plt.show()
