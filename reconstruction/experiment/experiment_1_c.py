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
from sklearn.neighbors import KDTree

from csv_modules.csv_writer import write_in_file
from experiment.utils_general import remove_get_dirs
from general_utils.database_utils import get_chains_pdb_db, get_graph_pdb_db, get_zd_pdb_db, get_zd_chain_pdb_db, \
  get_zd_chains_pdb_db
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_similar_pdb_struct, get_similar_pdb_chain_structural, get_ignore_pdbs, \
  get_similar_pdb_chain_sequential, get_percentage_pbs_check_file
from process_graph.graph_algorithm import graph_aligning
from process_mrc.miscellaneous import get_center_point_by_graph

headers = ['Pdb', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
           'Point Original syn dis', 'Point Test syn dis',
           'Resolution', 'Match',
           'Nodes father', 'Nodes test', 'Alignment note',
           'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']


def do_parallel_test(path_data,
                     result_cvs,
                     resolution_range=[5.0, 5.0], can_elements=None,
                     ignore_pdbs=[], percentage_data_set=10, file_checkpoint='check_expe_1c.pkl',
                     error_file='error.txt',
                     add_to_ignore_files=False,
                     can_groups=3,
                     min_can_chains=6):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      all_names = get_percentage_pbs_check_file(percentage_data_set, file_checkpoint, executor, min_can_chains)
      # all_names = ['6x0m']
      path = os.path.abspath(path_data)

      if not os.path.isdir(path):
        os.mkdir(path)

      complete_pdb = remove_get_dirs(path_data, can_csv=1, add_to_ignore_files=add_to_ignore_files)
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
        parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_aux, path, pdb_name,
                                                        result_cvs,
                                                        resolution, can_groups),

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


def do_parallel_test_aux(path, pdb_name, result_cvs, resolution, can_groups):
  local_path = os.path.join(path, pdb_name)
  if not os.path.exists(local_path):
    os.makedirs(local_path)
  os.makedirs(os.path.join(local_path, "flat"))

  # Ok with can_groups
  can_groups = max(1, can_groups)
  can_groups = min(len(get_chains_pdb_db(pdb_name)), can_groups)

  # Order points in tree
  list_chain = []
  original_graph = get_graph_pdb_db(pdb_name, resolution)
  data_kd = []
  dicc_node_pos = {}
  for i in original_graph.nodes:
    node_data = original_graph.nodes[i]
    x = node_data["cube_xyz_can"][0][0]
    y = node_data["cube_xyz_can"][0][1]
    z = node_data["cube_xyz_can"][0][2]
    total_div_points = node_data["cube_xyz_can"][1]

    x //= total_div_points
    y //= total_div_points
    z //= total_div_points

    list_chain.append(i)
    data_kd.append([x, y, z])
    dicc_node_pos[i] = [x, y, z]

  kd_tree = KDTree(np.array(data_kd))

  # Gen initial elements of groups
  max_distance_elements = []
  min_max_distance = 0
  for i in list_chain:
    dist, ind = kd_tree.query(np.array([dicc_node_pos[i]]), k=original_graph.number_of_nodes())
    dist = dist.tolist()
    ind = ind.tolist()
    dist = dist[0]
    ind = ind[0]

    if dist[(-can_groups) + 1] > min_max_distance:
      max_distance_elements = []
      min_max_distance = dist[(-can_groups) + 1]
      for k in range((-can_groups) + 1, 0, 1):
        max_distance_elements.append(list_chain[ind[k]])
      max_distance_elements.append(i)

  # Gen groups
  match_elements = max_distance_elements
  groups = np.expand_dims(max_distance_elements, axis=1).tolist()

  while len(match_elements) < original_graph.number_of_nodes():
    for i in range(can_groups):
      dist, ind = kd_tree.query(np.array([dicc_node_pos[groups[i][0]]]),
                                k=original_graph.number_of_nodes())
      dist = dist.tolist()
      ind = ind.tolist()
      dist = dist[0]
      ind = ind[0]
      for k in ind:
        if list_chain[ind[k]] not in match_elements:
          groups[i].append(list_chain[ind[k]])
          match_elements.append(list_chain[ind[k]])
          break

      if len(match_elements) >= original_graph.number_of_nodes():
        break

  # Do experiments
  for i in groups:
    father_graph = original_graph.copy()
    test_graph = original_graph.copy()

    nodes_to_delete = np.setdiff1d(list(father_graph.nodes), i).tolist()
    for k in nodes_to_delete:
      test_graph.remove_node(k)

    do_test(pdb_name, local_path, father_graph, test_graph, resolution, result_cvs)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.find('.') == -1 or directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      if os.path.isdir(path_remove):
        os.rmdir(path_remove)
      else:
        os.remove(path_remove)


def do_test(pdb_name, local_path, father_graph, test_graph, resolution, result_cvs):
  # Match syntetic
  original_match_list = [[i, i] for i in list(father_graph.nodes)]
  graph1_match_index = get_element_list(0, original_match_list)

  test_match_list = [[i, i] for i in list(test_graph.nodes)]
  graph2_match_index = get_element_list(1, test_match_list)

  start_time = time.time()
  center_point1 = get_center_point_by_graph(graph1_match_index, father_graph)
  center_point2 = get_center_point_by_graph(graph2_match_index, test_graph)
  time_center = time.time() - start_time

  # Real match process
  start_time = time.time()
  alignment_note, result = graph_aligning(father_graph, test_graph, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point_by_graph(graph1_match_index, father_graph)
    center_point2_1 = get_center_point_by_graph(graph2_match_index, test_graph)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[pdb_name, get_chains_pdb_db(pdb_name),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution, result,
                 list(father_graph.nodes), list(test_graph.nodes), alignment_note,
                 0,
                 time_center, time_center_test,
                 0,
                 time_aligning, 0]]

  write_in_file('{0}/{1}'.format(local_path, result_cvs), headers, data_write)
