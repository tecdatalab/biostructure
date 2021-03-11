import os
import random
import time
import traceback

import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from experiment.utils_experiment_1 import gen_keys_experiemnts
from experiment.utils_general import remove_get_dirs
from general_utils.database_utils import get_graph_pdb_db
from general_utils.download_utils import download_pdb
from general_utils.graph_utils import graph_is_connect
from general_utils.list_utils import get_element_list, combinations_i2jInK
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_ignore_pdbs, get_chains_pdb, get_all_pdb_name, get_percentage_pbs_check_file
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from process_mrc.miscellaneous import get_center_point, get_center_point_by_graph


def do_parallel_test_a(path_data, result_cvs_file, resolution_range=[5.0, 5.0], can_elements=None,
                       ignore_pdbs=[], range_incompleteness=[10.0, 15.0], percentage_data_set=10,
                       can_experiments_to_do=-1,
                       file_checkpoint='check_expe_1.pkl', add_to_ignore_files=False, error_file='error.txt'):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      # Get Pdbs to work
      all_names = get_percentage_pbs_check_file(percentage_data_set, file_checkpoint, executor)
      path = os.path.abspath(path_data)

      print("Before get pdbs work")

      if not os.path.isdir(path):
        os.mkdir(path)

      complete_pdb = remove_get_dirs(path_data, add_to_ignore_files=add_to_ignore_files)
      ignore_pdbs += complete_pdb

      # Add ignore files
      ignore_pdbs += get_ignore_pdbs()

      if can_elements is None:
        can_elements = len(all_names)

      parallel_jobs = []

      all_names = np.setdiff1d(np.array(all_names), np.array(ignore_pdbs)).tolist()[:can_elements]
      print("Do ", len(all_names), flush=True)
      for pdb_name in all_names:

        resolution = random.choices(resolution_range)[0]

        parallel_jobs.append([pdb_name,
                              executor.submit(do_parallel_test_a_aux, path, pdb_name, result_cvs_file,
                                              resolution, range_incompleteness, can_experiments_to_do),
                              resolution])
        # do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution, range_incompleteness, can_try_experiments)
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


def do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution, range_incompleteness,can_experiments_to_do):
  local_path = os.path.join(path, pdb_name)
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  graph_pdb = get_graph_pdb_db(pdb_name, resolution)

  headers_csv = ['Pdb', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                 'Point Original syn dis', 'Point Test syn dis',
                 'Resolution', 'Match note',
                 'Match', 'Father Chains', 'Test Chains',
                 'Time center', 'Time test', 'Time alignment',
                 'Graph original is connect', 'Graph test is connect']

  remove_try = combinations_i2jInK(graph_pdb.number_of_nodes(),
                                   round(graph_pdb.number_of_nodes() * (range_incompleteness[0]/100)),
                                   round(graph_pdb.number_of_nodes() * (range_incompleteness[1]/100)))
  graph_pdb_nodes = list(graph_pdb.nodes)
  random.shuffle(remove_try)

  for i in remove_try[:can_experiments_to_do]:

    graph_pdb_test = graph_pdb.copy()

    for node_pos in i:
      graph_pdb_test.remove_node(graph_pdb_nodes[node_pos])

    do_test(pdb_name, headers_csv, result_cvs_file, resolution, local_path,
            graph_pdb, graph_pdb_test)


  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_test(pdb_name, headers_csv, result_cvs_file, resolution, path_write,
              graph_pdb, graph_pdb_test):

  test_match_list = [[i, i] for i in list(graph_pdb_test.nodes)]

  graph1_match_index = get_element_list(0, test_match_list)
  graph2_match_index = get_element_list(1, test_match_list)

  start_time = time.time()
  center_point1 = get_center_point_by_graph(graph1_match_index, graph_pdb)
  center_point2 = get_center_point_by_graph(graph2_match_index, graph_pdb_test)
  time_center = time.time() - start_time

  start_time = time.time()
  alignment_note, result = graph_aligning(graph_pdb, graph_pdb_test, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point_by_graph(graph1_match_index, graph_pdb)
    center_point2_1 = get_center_point_by_graph(graph2_match_index, graph_pdb_test)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[pdb_name, list(graph_pdb.nodes),
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution, alignment_note,
                 result, list(graph_pdb.nodes), list(graph_pdb_test.nodes),
                 time_center, time_center_test,
                 time_aligning,
                 graph_is_connect(graph_pdb), graph_is_connect(graph_pdb_test)]]

  write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, data_write)
