import os
import pickle
import random
import shutil
import time
import traceback
import math

import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from sklearn.metrics import mean_squared_error

from csv_modules.csv_writer import write_in_file
from experiment.utils_general import remove_get_dirs, pdb_percentage
from general_utils.download_utils import download_pdb
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_similar_pdb_struct, get_similar_pdb_chain_structural, get_ignore_pdbs, \
  get_chains_pdb, get_similar_pdb_chain_sequential
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from process_mrc.miscellaneous import get_center_point

headers_chain = ['Pdb', 'Pdb work', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                 'Point Original syn dis', 'Point Test syn dis',
                 'Resolution', 'Match',
                 'Score PDB', '3D ZD note', 'Alignment note',
                 'Changed chain', 'Changed chain Pdb work',
                 'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

headers_struct = ['Pdb', 'Pdb work', 'Chains', 'Work Chains', 'Point Original', 'Point Test', 'Point Original syn',
                  'Point Test syn',
                  'Point Original syn dis', 'Point Test syn dis',
                  'Resolution', 'Match',
                  'Score PDB', '3D ZD note (complete)', 'Alignment note',
                  'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']


def do_parallel_test(path_data,
                     result_cvs_chain, result_cvs_struct, result_cvs_secuencial,
                     resolution_range=[5.0, 5.0], can_elements=None,
                     ignore_pdbs=[], percentage_data_set=10, file_checkpoint='check_expe_1a.pkl',
                     error_file='error.txt',
                     can_chain_test=3, can_struct_test=3, can_secuencial_test=3,
                     add_to_ignore_files=False):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      # if not os.path.exists(file_checkpoint):
      #   all_names = pdb_percentage(percentage_data_set, executor)  # 169315
      #   open_file = open(file_checkpoint, "wb")
      #   pickle.dump(all_names, open_file)
      #   open_file.close()
      # else:
      #   open_file = open(file_checkpoint, "rb")
      #   all_names = pickle.load(open_file)
      #   open_file.close()

      # print(all_names, flush=True)
      # all_names = ['1aig']
      # all_names = ['1q8l']
      # all_names = ['1a0t']
      all_names = ['109l']
      print("Before get pdb names")

      path = os.path.abspath(path_data)
      # print(path)

      if not os.path.isdir(path):
        os.mkdir(path)

      complete_pdb = remove_get_dirs(path_data, can_csv=3, add_to_ignore_files=add_to_ignore_files)
      ignore_pdbs += complete_pdb
      # Add ignore files
      ignore_pdbs += get_ignore_pdbs()

      if can_elements is None:
        can_elements = len(all_names)

      parallel_jobs = []

      all_names = np.setdiff1d(np.array(all_names), np.array(ignore_pdbs)).tolist()[:can_elements]
      print("Do ", len(all_names), flush=True)
      for pdb_name in all_names:
        resolution = random.uniform(resolution_range[0], resolution_range[1])
        # resolution = 3.8680

        # print(pdb_name, con2/can_elements)
        parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_aux, path, pdb_name, result_cvs_chain,
                                                        result_cvs_struct, result_cvs_secuencial,
                                                        resolution,
                                                        can_chain_test, can_struct_test, can_secuencial_test),

                              resolution])
        # do_parallel_test_a_aux(path, pdb_name, result_cvs_chain, result_cvs_struct, resolution, can_chain_test,
        # can_struct_test)
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          with open(error_file, "a+") as myfile:
            myfile.write(f[0])
            myfile.write("\n")
            myfile.write(str(f[2]))
            myfile.write("\n")
            myfile.write(str(type(e).__name__))
            myfile.write("\n")
            myfile.write(str(e))
            myfile.write("\n")
            myfile.write(str(traceback.format_exc()))
            myfile.write("\n\n\n\n")


def do_parallel_test_aux(path, pdb_name, result_cvs_chain, result_cvs_struct, result_cvs_secuencial, resolution,
                         can_chain_test, can_struct_test, can_secuencial_test):
  local_path = path + "/" + pdb_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)
  path_of_pdb = '{0}/{1}.pdb'.format(local_path, pdb_name)

  download_pdb(pdb_name, path_of_pdb)
  chains = get_chains_pdb(path_of_pdb)

  result_struct = []
  temp = get_similar_pdb_struct(pdb_name, -1)
  for i in temp:
    add_data = [pdb_name, i[0], i[1]]
    result_struct.append(add_data)

  result_chain_struct = []
  for chain in chains:
    temp = get_similar_pdb_chain_structural(pdb_name, chain, -1)
    for i in temp:
      add_data = [pdb_name, i[0], i[1], chain]
      result_chain_struct.append(add_data)

  result_chain_sequence = []
  for chain in chains:
    temp = get_similar_pdb_chain_sequential(pdb_name, chain, -1)
    for i in temp:
      add_data = [pdb_name, i[0], i[1], chain]
      result_chain_sequence.append(add_data)

  # Clen not do
  result_struct = get_experiments_to_do(result_struct, can_struct_test)
  result_chain_struct = get_experiments_to_do(result_chain_struct, can_chain_test)
  result_chain_sequence = get_experiments_to_do(result_chain_sequence, can_secuencial_test)

  for struct_score in result_struct:
    do_test_struct(local_path, struct_score, resolution, result_cvs_struct)
  if result_struct==[]:
    do_test_struct(local_path, None, resolution, result_cvs_struct)

  for chain_score in result_chain_struct:
    do_test_chain(local_path, chain_score, resolution, result_cvs_chain)
  if result_chain_struct==[]:
    do_test_chain(local_path, None, resolution, result_cvs_chain)

  for sequential_score in result_chain_sequence:
    do_test_chain(local_path, sequential_score, resolution, result_cvs_secuencial)
  if result_chain_sequence == []:
    do_test_chain(local_path, None, resolution, result_cvs_secuencial)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def get_experiments_to_do(list_possibles, cant_by_range):
  dic_values = {}
  for exp in list_possibles:
    key = math.floor(exp[2] * 10)
    if dic_values.get(key) == None:
      dic_values[key] = [exp]
    else:
      dic_values[key].append(exp)

  result = []

  for key in dic_values.keys():
    random.choice(dic_values[key])
    result += dic_values[key][:cant_by_range]
  return result


def do_test_struct(path, struct_score, resolution, result_cvs_struct):
  if struct_score == None:
    write_in_file('{0}/{1}'.format(path, result_cvs_struct), headers_struct, [[]])
    return

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
                 resolution, result,
                 score_pdb_work, mean_squared_error(pdb_work_mrc.zd_descriptors, pdb_test_mrc.zd_descriptors),
                 alignment_note,
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(path, result_cvs_struct), headers_struct, data_write)
  shutil.rmtree(local_path)


def do_test_chain(path, chain_score, resolution, result_cvs_chain):
  if chain_score == None:
    write_in_file('{0}/{1}'.format(path, result_cvs_chain), headers_chain, [[]])
    return

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

  pdb_work, score_chain_work, chain_work = chain_score[1], chain_score[2], chain_score[3]
  id_work = chains_to_segment[chain_work].id_segment

  # Change for real method
  zernike_zd_descriptors_chain_pdb, score_chain_changed, chain_changed = get_zd_descriptors_chain_pdb(pdb_work,
                                                                                                      chains_to_segment[
                                                                                                        chain_work].zd_descriptors,
                                                                                                      local_path,
                                                                                                      resolution)

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
                 resolution, result,
                 score_chain_work, score_chain_changed, alignment_note,
                 chain_work, chain_changed,
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(path, result_cvs_chain), headers_chain, data_write)
  shutil.rmtree(local_path)


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
