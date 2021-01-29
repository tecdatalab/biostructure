import os
import pickle
import random
import shutil
import time
import copy
import numpy as np
from sklearn.metrics import mean_squared_error
import traceback

from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from experiment.utils_experiment_1 import gen_keys_experiemnts
from experiment.utils_general import remove_get_dirs, pdb_percentage
from general_utils.download_utils import download_pdb
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_similar_pdb_struct, get_similar_pdb_chain, get_pdb_no_work, get_ignore_pdbs, \
    get_chains_pdb, get_all_pdb_name
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from process_mrc.miscellaneous import get_center_point


def do_parallel_test_a(path_data, result_cvs_chain, result_cvs_struct, result_cvs_secuencial,
                       all_name_struct, all_name_chain, all_name_secuencial,
                       resolution_range=[5.0, 5.0]):

  path = os.path.abspath(path_data)
  if not os.path.isdir(path):
    os.mkdir(path)

  for struct_score in all_name_struct:
    resolution = random.uniform(resolution_range[0], resolution_range[1])
    do_test_struct(path_data, struct_score, resolution, result_cvs_struct)

  for chain_score in all_name_chain:
    resolution = random.uniform(resolution_range[0], resolution_range[1])
    do_test_chain(path_data, chain_score, resolution, result_cvs_chain)

  for secuencial_score in all_name_secuencial:
    resolution = random.uniform(resolution_range[0], resolution_range[1])
    do_test_secuencial(path_data, secuencial_score, resolution, result_cvs_secuencial)


def do_test_struct(path, struct_score, resolution, result_cvs_struct):
  local_path = path + "/" + struct_score[0]
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(struct_score[0], '{0}/{1}.pdb'.format(local_path, struct_score[0]))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(local_path, struct_score[0]))

  start_time = time.time()
  pdb_to_mrc_chains(True, False, resolution, '{0}/{1}.pdb'.format(local_path, struct_score[0]), path, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, struct_score[0]))
  time_eman = time.time() - start_time

  chains_to_segment = {}
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, struct_score[0], chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    chains_to_segment[chain] = segment
    con_id_segment += 1
  time_segment = time.time() - start_time

  original_segments = []
  for i in chains_to_segment.keys():
    original_segments.append(chains_to_segment[i])

  headers_struct = ['Pdb', 'Pdb work', 'Chains', 'Work Chains', 'Point Original', 'Point Test', 'Point Original syn',
                    'Point Test syn',
                    'Point Original syn dis', 'Point Test syn dis',
                    'Resolution',
                    'Score pdb work', 'ZD Complete note','Alinament note',
                    'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  # Change for real method
  pdb_work, score_pdb_work, test_segments, work_chains = get_segments_struct_test(struct_score[1:], local_path,
                                                                                  resolution)

  original_match_list = [[i.id_segment, i.id_segment] for i in original_segments]
  test_match_list = [[i.id_segment, i.id_segment] for i in test_segments]

  graph1_match_index = get_element_list(0, original_match_list)
  graph2_match_index = get_element_list(1, test_match_list)

  start_time = time.time()
  center_point1 = get_center_point(graph1_match_index, original_segments, 0)
  center_point2 = get_center_point(graph2_match_index, test_segments, 0)
  time_center = time.time() - start_time

  # Generate test syntetic
  start_time = time.time()
  graph1 = generate_graph(original_segments, 50, 0, 6, 1)
  graph2 = generate_graph(test_segments, 50, 0, 6, 1)
  time_graph = time.time() - start_time

  start_time = time.time()
  alignment_note, result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point(graph1_match_index, original_segments, 0)
    center_point2_1 = get_center_point(graph2_match_index, test_segments, 0)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  download_pdb(struct_score[1], '{0}/{1}.pdb'.format(local_path, struct_score[1]))
  pdb_to_mrc_chains(True, False, resolution, '{0}/{1}.pdb'.format(local_path, struct_score[1]), local_path)

  pdb_work_mrc = get_mrc_one('{0}/{1}.mrc'.format(local_path, struct_score[0]))[0][0]
  pdb_test_mrc = get_mrc_one('{0}/{1}/{1}.mrc'.format(local_path, struct_score[1]))[0][0]
  shutil.rmtree(local_path + '/' + struct_score[1])

  data_write = [[struct_score[0], pdb_work, list(chains_to_segment.keys()), work_chains,
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution,
                 score_pdb_work, mean_squared_error(pdb_work_mrc.zd_descriptors, pdb_test_mrc.zd_descriptors),
                 alignment_note,
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(local_path, result_cvs_struct), headers_struct, data_write)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_test_chain(path, chain_score, resolution, result_cvs_chain):
  local_path = path + "/" + chain_score[0]
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(chain_score[0], '{0}/{1}.pdb'.format(local_path, chain_score[0]))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(local_path, chain_score[0]))

  start_time = time.time()
  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, chain_score[0]), path, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, chain_score[0]))
  time_eman = time.time() - start_time

  chains_to_segment = {}
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, chain_score[0], chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    chains_to_segment[chain] = segment
    con_id_segment += 1
  time_segment = time.time() - start_time

  original_segments = []
  for i in chains_to_segment.keys():
    original_segments.append(chains_to_segment[i])

  headers_chain = ['Pdb', 'Pdb work', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                   'Point Original syn dis', 'Point Test syn dis',
                   'Resolution',
                   'Score pdb work', 'Score chain changed', 'Alinament note',
                   'Changed chain', 'Changed chain Pdb work',
                   'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  pdb_work, score_chain_work, chain_work = chain_score[1], chain_score[2], chain_score[3]
  id_work = chains_to_segment[chain_work].id_segment

  # Change for real method
  zernike_zd_descriptors_chain_pdb, score_chain_changed, chain_changed = get_zd_descriptors_chain_pdb(pdb_work,
                                                                                       chains_to_segment[
                                                                                         chain_work].zd_descriptors,
                                                                                       local_path, resolution)

  original_match_list = [[i.id_segment, i.id_segment] for i in original_segments]

  graph1_match_index = get_element_list(0, original_match_list)
  graph2_match_index = get_element_list(1, original_match_list)

  start_time = time.time()
  center_point1 = get_center_point(graph1_match_index, original_segments, 0)
  center_point2 = get_center_point(graph2_match_index, original_segments, 0)
  time_center = time.time() - start_time

  # Generate test syntetic
  start_time = time.time()
  graph1 = generate_graph(original_segments, 50, 0, 6, 1)
  graph2 = graph1.copy()
  graph2.nodes[id_work]["zd_descriptors"] = zernike_zd_descriptors_chain_pdb
  time_graph = time.time() - start_time

  start_time = time.time()
  alignment_note, result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point(graph1_match_index, original_segments, 0)
    center_point2_1 = get_center_point(graph2_match_index, original_segments, 0)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[chain_score[0], pdb_work, list(chains_to_segment.keys()),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution,
                 score_chain_work, score_chain_changed, alignment_note,
                 chain_work, chain_changed,
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(local_path, result_cvs_chain), headers_chain, data_write)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_test_secuencial(path, secuencial_score, resolution, result_cvs_chain):
  local_path = path + "/" + secuencial_score[0]
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(secuencial_score[0], '{0}/{1}.pdb'.format(local_path, secuencial_score[0]))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(local_path, secuencial_score[0]))

  start_time = time.time()
  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, secuencial_score[0]), path, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, secuencial_score[0]))
  time_eman = time.time() - start_time

  chains_to_segment = {}
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, secuencial_score[0], chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    chains_to_segment[chain] = segment
    con_id_segment += 1
  time_segment = time.time() - start_time

  original_segments = []
  for i in chains_to_segment.keys():
    original_segments.append(chains_to_segment[i])

  headers_chain = ['Pdb', 'Pdb work', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                   'Point Original syn dis', 'Point Test syn dis',
                   'Resolution',
                   'Score pdb work', 'Score chain changed', 'Alinament note',
                   'Changed chain', 'Changed chain Pdb work',
                   'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  pdb_work, score_chain_work, chain_work = secuencial_score[1], secuencial_score[2], secuencial_score[3]
  id_work = chains_to_segment[chain_work].id_segment

  # Change for real method
  zernike_zd_descriptors_chain_pdb, score_chain_changed = get_zd_descriptors_one_chain_pdb(pdb_work,
                                                                                       chains_to_segment[
                                                                                         chain_work].zd_descriptors,
                                                                                       local_path, resolution,
                                                                                           secuencial_score[4])

  original_match_list = [[i.id_segment, i.id_segment] for i in original_segments]

  graph1_match_index = get_element_list(0, original_match_list)
  graph2_match_index = get_element_list(1, original_match_list)

  start_time = time.time()
  center_point1 = get_center_point(graph1_match_index, original_segments, 0)
  center_point2 = get_center_point(graph2_match_index, original_segments, 0)
  time_center = time.time() - start_time

  # Generate test syntetic
  start_time = time.time()
  graph1 = generate_graph(original_segments, 50, 0, 6, 1)
  graph2 = graph1.copy()
  graph2.nodes[id_work]["zd_descriptors"] = zernike_zd_descriptors_chain_pdb
  time_graph = time.time() - start_time

  start_time = time.time()
  alignment_note, result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point(graph1_match_index, original_segments, 0)
    center_point2_1 = get_center_point(graph2_match_index, original_segments, 0)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[secuencial_score[0], pdb_work, list(chains_to_segment.keys()),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution,
                 score_chain_work, score_chain_changed, alignment_note,
                 chain_work, secuencial_score[4],
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(local_path, result_cvs_chain), headers_chain, data_write)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)

def get_segments_struct_test(pdb_score, path_work, resolution):
  download_pdb(pdb_score[0], '{0}/{1}.pdb'.format(path_work, pdb_score[0]))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_work, pdb_score[0]))

  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(path_work, pdb_score[0]), path_work, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(path_work, pdb_score[0]))

  all_segments = []
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one(
      '{0}/{1}_{2}.mrc'.format(path_work + '/' + pdb_score[0], pdb_score[0], chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    all_segments.append(segment)
    con_id_segment += 1

  shutil.rmtree(path_work + '/' + pdb_score[0])

  return pdb_score[0], pdb_score[1], all_segments, chains


def get_zd_descriptors_chain_pdb(pdb_work, zd_descriptors_compare, path_work, resolution):
  download_pdb(pdb_work, '{0}/{1}.pdb'.format(path_work, pdb_work))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_work, pdb_work))

  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(path_work, pdb_work), path_work, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(path_work, pdb_work))

  best_distance = float('inf')
  result = []
  chain_changed = ''

  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one(
      '{0}/{1}_{2}.mrc'.format(path_work + '/' + pdb_work, pdb_work, chain))

    mse = mean_squared_error(zd_descriptors_compare, segments_graph_simulate[0].zd_descriptors)
    if mse < best_distance:
      best_distance = mse
      result = segments_graph_simulate[0].zd_descriptors
      chain_changed = chain

  shutil.rmtree(path_work + '/' + pdb_work)

  return result, best_distance, chain_changed


def get_zd_descriptors_one_chain_pdb(pdb_work, zd_descriptors_compare, path_work, resolution, chain):
  download_pdb(pdb_work, '{0}/{1}.pdb'.format(path_work, pdb_work))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_work, pdb_work))

  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(path_work, pdb_work), path_work, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(path_work, pdb_work))



  segments_graph_simulate, _ = get_mrc_one(
      '{0}/{1}_{2}.mrc'.format(path_work + '/' + pdb_work, pdb_work, chain))

  best_distance = mean_squared_error(zd_descriptors_compare, segments_graph_simulate[0].zd_descriptors)
  result = segments_graph_simulate[0].zd_descriptors

  shutil.rmtree(path_work + '/' + pdb_work)

  return result, best_distance


def experiment_1_a():
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"

  all_name_struct = [
    ["4udt", "2pyf", "0.813041882029752"],
    ["1ku8", "4r6n", "0.4379816205362538"],
    ["1ku8", "6kqb", "0.13073256037178124"]
  ]
  all_name_chain = [
    ["3fns", "2igx", "0.36695932635646", "B"],
    ["1kv3", "6kzb", "0.38971299062784465", "A"],
    ["1kv3", "2ygb", "0.09056548598533672", "A"]
  ]
  #La primera cadena es la propia
  #La segunda cadena es la cambiada
  all_name_secuencial = [
    ["1mbn", "2mye", str(float('1.31E+01')), "A", "A"],
    ["1kxq", "1kxt", str(float('100')), "A", "E"],
    ["1kxq", "3bai", str(float('86.9')), "A", "A"],
    ["1kxq", "1kcl", str(float('23.5')), "A", "A"],

  ]
  # local_path = "/work/lcastillo"
  print("Start force")
  do_parallel_test_a("{0}/data_experiment_1_a_a_v1".format(local_path),
                     result_cvs_chain="result_chain.csv",
                     result_cvs_struct="result_struct.csv",
                     result_cvs_secuencial="result_secuencial.csv",
                     all_name_struct=all_name_struct,
                     all_name_chain=all_name_chain,
                     all_name_secuencial=all_name_secuencial,
                     resolution_range=[3.5, 9.5])

  print("Finish")


if __name__ == '__main__':
  experiment_1_a()
  #result =get_similar_pdb_struct("1KU8", 100)
  #result = get_similar_pdb_chain("1kv3", "A", 100)
  #print(result)
