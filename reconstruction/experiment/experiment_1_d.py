import os
import random
import time
import traceback
import tempfile
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from experiment.utils_general import remove_get_dirs
from general_utils.database_utils import get_graph_pdb_db, get_chains_pdb_db
from general_utils.graph_utils import graph_is_connect
from general_utils.list_utils import get_element_list, combinations_i2jInK
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_ignore_pdbs, get_percentage_pbs_check_file
from process_graph.graph_algorithm import graph_aligning
from process_mrc.miscellaneous import get_center_point_by_graph


def do_parallel_test_a(path_data, result_cvs_file, resolution_range=[5.0, 5.0], pdbs_work=[],
                       error_file='error.txt'):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      path = os.path.abspath(path_data)

      if not os.path.isdir(path):
        os.mkdir(path)

      parallel_jobs = []

      for pdb_name in pdbs_work:
        actual_resolution = random.choice(resolution_range)
        total_combinations = combinations_i2jInK(len(get_chains_pdb_db(pdb_name)), 1,
                                                 len(get_chains_pdb_db(pdb_name))-1, check_max=False)

        for combination_to_do_delete in total_combinations:
          parallel_jobs.append([pdb_name,
                                executor.submit(do_parallel_test_aux, path, pdb_name, result_cvs_file,
                                                actual_resolution, combination_to_do_delete),
                                actual_resolution])
        # do_parallel_test_a_aux
      total_experiments_do = len(parallel_jobs)
      con = 0
      for f in parallel_jobs:
        con+=1
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
        print("Actual:", (con/total_experiments_do)*100, "%", flush=True)



def do_parallel_test_aux(path, pdb_name, result_cvs_file, resolution, combination_to_do_delete):
  local_path = tempfile.mkdtemp(dir=path)

  graph_pdb = get_graph_pdb_db(pdb_name, resolution)

  headers_csv = ['Pdb', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                 'Point Original syn dis', 'Point Test syn dis',
                 'Resolution', 'Match note',
                 'Match', 'Father Chains', 'Test Chains',
                 'Time center', 'Time test', 'Time alignment',
                 'Graph original is connect', 'Graph test is connect']

  graph_pdb_nodes = list(graph_pdb.nodes)
  graph_pdb_test = graph_pdb.copy()

  for node_pos in combination_to_do_delete:
    graph_pdb_test.remove_node(graph_pdb_nodes[node_pos])

  do_test(pdb_name, headers_csv, result_cvs_file, resolution, local_path,
          graph_pdb, graph_pdb_test)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.find('.') == -1 or directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      if os.path.isdir(path_remove):
        os.rmdir(path_remove)
      else:
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
