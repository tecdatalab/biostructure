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
from experiment.utils_general import remove_get_dirs, check_RMSD_result_all, check_RMSD_result_algorithm
from general_utils.database_utils import get_chains_pdb_db, get_graph_pdb_db, get_zd_pdb_db, get_zd_chain_pdb_db, \
  get_zd_chains_pdb_db
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_similar_pdb_struct, get_similar_pdb_chain_structural, get_ignore_pdbs, \
  get_similar_pdb_chain_sequential, get_percentage_pbs_check_file
from general_utils.temp_utils import gen_dir, free_dir
from process_graph.graph_algorithm import graph_aligning
from process_mrc.miscellaneous import get_center_point_by_graph

headers_chain = ['Pdb', 'Pdb work', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                 'Point Original syn dis', 'Point Test syn dis',
                 'Resolution', 'Match',
                 'Score PDB', '3D ZD note', 'Alignment note',
                 'Changed chain', 'Changed chain Pdb work',
                 'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2',
                 'RMSD']

headers_struct = ['Pdb', 'Pdb work', 'Chains', 'Work Chains', 'Point Original', 'Point Test', 'Point Original syn',
                  'Point Test syn',
                  'Point Original syn dis', 'Point Test syn dis',
                  'Resolution', 'Match',
                  'Score PDB', '3D ZD note (complete)', 'Alignment note',
                  'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2',
                  'RMSD']


def do_parallel_test(path_data,
                     result_cvs_chain, result_cvs_struct, result_cvs_secuencial,
                     resolution_range=[5.0, 5.0], can_elements=None,
                     ignore_pdbs=[], percentage_data_set=10, file_checkpoint='check_expe_1b.pkl',
                     error_file='error_log_expe_1b.txt',
                     can_chain_test=3, can_struct_test=3, can_secuencial_test=3,
                     add_to_ignore_files=False):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      all_names = get_percentage_pbs_check_file(percentage_data_set, file_checkpoint, executor)
      # all_names = ['4u0d']
      path = os.path.abspath(path_data)

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
      random.shuffle(all_names)
      print("Do ", len(all_names), flush=True)
      for pdb_name in all_names:
        resolution = random.choices(resolution_range)[0]
        # resolution = 10

        # print(pdb_name, con2/can_elements)
        parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_aux, path, pdb_name, result_cvs_chain,
                                                        result_cvs_struct, result_cvs_secuencial,
                                                        resolution,
                                                        can_chain_test, can_struct_test, can_secuencial_test,
                                                        error_file),

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
                         can_chain_test, can_struct_test, can_secuencial_test, error_file):
  try:
    local_path = os.path.join(path, pdb_name)
    if not os.path.exists(local_path):
      os.makedirs(local_path)

    RMSD_dir = gen_dir()
    chains = get_chains_pdb_db(pdb_name)

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
      do_test_struct(local_path, struct_score, resolution, result_cvs_struct, RMSD_dir)
    if result_struct == []:
      do_test_struct(local_path, None, resolution, result_cvs_struct, RMSD_dir)

    for chain_score in result_chain_struct:
      do_test_chain(local_path, chain_score, resolution, result_cvs_chain, RMSD_dir)
    if result_chain_struct == []:
      do_test_chain(local_path, None, resolution, result_cvs_chain, RMSD_dir)

    for sequential_score in result_chain_sequence:
      do_test_chain(local_path, sequential_score, resolution, result_cvs_secuencial, RMSD_dir)
    if result_chain_sequence == []:
      do_test_chain(local_path, None, resolution, result_cvs_secuencial, RMSD_dir)

    dirs = os.listdir(local_path)
    free_dir(RMSD_dir)

    for directory in dirs:
      if directory.find('.') == -1 or directory.split('.')[1] != 'csv':
        path_remove = '{0}/{1}'.format(local_path, directory)
        if os.path.isdir(path_remove):
          os.rmdir(path_remove)
        else:
          os.remove(path_remove)
  except Exception as e:
    with open(error_file, "a+") as myfile:
      myfile.write("Not control error")
      myfile.write("\n")
      myfile.write(pdb_name)
      myfile.write("\n")
      myfile.write(str(resolution))
      myfile.write("\n")
      myfile.write(str(type(e).__name__))
      myfile.write("\n")
      myfile.write(str(e))
      myfile.write("\n")
      myfile.write(str(traceback.format_exc()))
      myfile.write("\n\n\n\n")


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
    random.shuffle(dic_values[key])
    result += dic_values[key][:cant_by_range]
  return result


def do_test_struct(path, struct_score, resolution, result_cvs_struct, RMSD_dir):
  if struct_score == None:
    write_in_file('{0}/{1}'.format(path, result_cvs_struct), headers_struct, [[]])
    return

  local_path = path + "/" + struct_score[0]
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  # Varaibles
  pdb_test = struct_score[0]
  pdb_work = struct_score[1]
  score_pdb_work = struct_score[2]

  # Generate graphs
  start_time = time.time()
  graph1 = get_graph_pdb_db(pdb_test, resolution)
  graph2 = get_graph_pdb_db(pdb_work, resolution)
  time_graph = time.time() - start_time

  # Match syntetic
  original_match_list = [[i, i] for i in list(graph1.nodes)]
  test_match_list = [[i, i] for i in list(graph2.nodes)]

  graph1_match_index = get_element_list(0, original_match_list)
  graph2_match_index = get_element_list(1, test_match_list)

  start_time = time.time()
  center_point1 = get_center_point_by_graph(graph1_match_index, graph1)
  center_point2 = get_center_point_by_graph(graph2_match_index, graph2)
  time_center = time.time() - start_time

  # Real match process
  start_time = time.time()
  alignment_note, result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point_by_graph(graph1_match_index, graph1)
    center_point2_1 = get_center_point_by_graph(graph2_match_index, graph2)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[pdb_test, pdb_work, get_chains_pdb_db(pdb_test), get_chains_pdb_db(pdb_work),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution, result,
                 score_pdb_work, mean_squared_error(get_zd_pdb_db(pdb_work, resolution),
                                                    get_zd_pdb_db(pdb_test, resolution)),
                 alignment_note,
                 0,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, 0,
                 check_RMSD_result_all(RMSD_dir,
                                       get_chains_pdb_db(pdb_test),
                                       get_chains_pdb_db(pdb_work),
                                       result,
                                       pdb_test,
                                       pdb_work)
                 ]]

  write_in_file('{0}/{1}'.format(path, result_cvs_struct), headers_struct, data_write)
  shutil.rmtree(local_path)


def do_test_chain(path, chain_score, resolution, result_cvs_chain, RMSD_dir):
  if chain_score == None:
    write_in_file('{0}/{1}'.format(path, result_cvs_chain), headers_chain, [[]])
    return

  local_path = os.path.join(path, chain_score[0])
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  # Varaibles
  pdb_test = chain_score[0]
  pdb_work = chain_score[1]
  score_chain_work = chain_score[2]
  chain_work = chain_score[3]

  # Generate First graph (original)
  graph1 = get_graph_pdb_db(pdb_test, resolution)

  # Get pos element in graph
  id_work = get_chains_pdb_db(pdb_test).index(chain_work) + 1

  # Change for real method
  zernike_zd_descriptors_chain_pdb, score_chain_changed, chain_changed = \
    get_zd_descriptors_chain_pdb(pdb_work,
                                 get_zd_chain_pdb_db(pdb_test, chain_work, resolution),
                                 resolution)

  # Match syntetic
  original_match_list = [[i, i] for i in list(graph1.nodes)]
  graph1_match_index = get_element_list(0, original_match_list)
  graph2_match_index = get_element_list(1, original_match_list)

  start_time = time.time()
  center_point1 = get_center_point_by_graph(graph1_match_index, graph1)
  center_point2 = get_center_point_by_graph(graph2_match_index, graph1)
  time_center = time.time() - start_time

  # Change matcheg chain
  start_time = time.time()
  graph2 = graph1.copy()
  graph2.nodes[id_work]["zd_descriptors"] = zernike_zd_descriptors_chain_pdb
  time_graph = time.time() - start_time

  # Real match process
  start_time = time.time()
  alignment_note, result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point_by_graph(graph1_match_index, graph1)
    center_point2_1 = get_center_point_by_graph(graph2_match_index, graph2)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[chain_score[0], pdb_work, get_chains_pdb_db(pdb_test),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution, result,
                 score_chain_work, score_chain_changed, alignment_note,
                 chain_work, chain_changed,
                 0,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, 0,
                 check_RMSD_result_algorithm(RMSD_dir,
                                             get_chains_pdb_db(pdb_test),
                                             result,
                                             pdb_test,
                                             pdb_work,
                                             chain_work,
                                             chain_changed)]]

  write_in_file('{0}/{1}'.format(path, result_cvs_chain), headers_chain, data_write)
  shutil.rmtree(local_path)


def get_zd_descriptors_chain_pdb(pdb_work, zd_descriptors_compare, resolution):
  best_distance = float('inf')
  result = None
  chain_changed = ''

  for chain_zd in get_zd_chains_pdb_db(pdb_work, resolution):
    mse = mean_squared_error(zd_descriptors_compare, chain_zd[1])
    if mse < best_distance:
      best_distance = mse
      result = chain_zd[1]
      chain_changed = chain_zd[0]

  return result, best_distance, chain_changed
