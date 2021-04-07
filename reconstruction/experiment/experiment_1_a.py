import os
import random
import time
import traceback

import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from sklearn.metrics import mean_squared_error

from csv_modules.csv_writer import write_in_file
from experiment.utils_general import remove_get_dirs
from general_utils.database_utils import get_chains_pdb_db, get_graph_pdb_db, get_zd_pdb_db
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from general_utils.pdb_utils import get_similar_pdb_struct, get_similar_pdb_chain_structural, get_ignore_pdbs, \
  get_percentage_pbs_check_file, get_similar_pdb_chain_sequential
from process_graph.graph_algorithm import graph_aligning
from process_mrc.miscellaneous import get_center_point_by_graph

header_csv = ['Pdb', 'Pdb work', 'Chains', 'Work Chains', 'Point Original', 'Point Test', 'Point Original syn',
              'Point Test syn',
              'Point Original syn dis', 'Point Test syn dis',
              'Resolution',
              'Score pdb work', 'Score ZD',
              'Time center', 'Time test', 'Time alignment']


def do_parallel_test_a(path_data, result_cvs_chain, result_cvs_sequence, result_cvs_struct,
                       resolution_range=[4, 4], can_elements=None, ignore_pdbs=[], percentage_data_set=10,
                       file_checkpoint='check_expe_1a.pkl', error_file='error.txt',
                       can_chain_test=3, can_sequence_test=3, can_struct_test=3,
                       add_to_ignore_files=False):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      # Get Pdbs to work
      all_names = get_percentage_pbs_check_file(percentage_data_set, file_checkpoint, executor)
      # all_names = ['1c5f']

      print("Before get pdb names")

      # Check path work
      path = os.path.abspath(path_data)
      if not os.path.isdir(path):
        os.mkdir(path)

      # Add ignore files
      complete_pdb = remove_get_dirs(path_data, can_csv=3, add_to_ignore_files=add_to_ignore_files)
      ignore_pdbs += complete_pdb
      ignore_pdbs += get_ignore_pdbs()

      if can_elements is None:
        can_elements = len(all_names)

      parallel_jobs = []

      # Remove ignore files in work files
      all_names = np.setdiff1d(np.array(all_names), np.array(ignore_pdbs)).tolist()[:can_elements]

      print("Do ", len(all_names), flush=True)
      for pdb_name in all_names:
        resolution = random.choices(resolution_range)[0]

        parallel_jobs.append([pdb_name,
                              executor.submit(do_parallel_test_a_aux, path, pdb_name,
                                              result_cvs_chain, result_cvs_sequence, result_cvs_struct,
                                              resolution,
                                              can_chain_test, can_sequence_test, can_struct_test),
                              resolution])
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


def do_parallel_test_a_aux(path, pdb_name, result_cvs_chain, result_cvs_sequence, result_cvs_struct, resolution,
                           can_chain_test, can_sequence_test, can_struct_test):
  local_path = os.path.join(path, pdb_name)
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  pdb_origin_graph = get_graph_pdb_db(pdb_name, resolution)
  chains = get_chains_pdb_db(pdb_name)

  # Get PDBs to do in experiment
  list_possibles_pdb_chain_struct = []
  for chain in chains:
    print(chain)
    temp = get_similar_pdb_chain_structural(pdb_name, chain, -1)
    list_possibles_pdb_chain_struct += temp


  list_possibles_pdb_chain_sequence = []
  for chain in chains:
    print(chain)
    temp = get_similar_pdb_chain_sequential(pdb_name, chain, -1)
    list_possibles_pdb_chain_sequence += temp

  list_possibles_pdb_struct = get_similar_pdb_struct(pdb_name, -1)

  # Do experiments

  do_test_pdb_list(pdb_name, header_csv, result_cvs_sequence, pdb_origin_graph, resolution, local_path,
                   can_sequence_test, list_possibles_pdb_chain_sequence)

  do_test_pdb_list(pdb_name, header_csv, result_cvs_chain, pdb_origin_graph, resolution, local_path,
                   can_chain_test, list_possibles_pdb_chain_struct)

  do_test_pdb_list(pdb_name, header_csv, result_cvs_struct, pdb_origin_graph, resolution, local_path,
                   can_struct_test, list_possibles_pdb_struct)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.find('.') == -1 or directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      if os.path.isdir(path_remove):
        os.rmdir(path_remove)
      else:
        os.remove(path_remove)


def do_test_pdb_list(pdb_name, headers_csv, result_cvs_file, pdb_origin_graph, resolution, path_write, can_chain_test,
                     list_possibles_pdb):
  random.shuffle(list_possibles_pdb)

  for i in range(min(can_chain_test, len(list_possibles_pdb))):

    # Change for real method
    pdb_work = list_possibles_pdb[i][0]
    score_pdb_work = list_possibles_pdb[i][1]
    graph_pdb_work = get_graph_pdb_db(pdb_work, resolution)

    original_match_list = [[i, i] for i in list(pdb_origin_graph.nodes)]
    test_match_list = [[i, i] for i in list(graph_pdb_work.nodes)]

    graph1_match_index = get_element_list(0, original_match_list)
    graph2_match_index = get_element_list(1, test_match_list)

    start_time = time.time()
    center_point1 = get_center_point_by_graph(graph1_match_index, pdb_origin_graph)
    center_point2 = get_center_point_by_graph(graph2_match_index, graph_pdb_work)
    time_center = time.time() - start_time

    # Generate test syntetic
    graph1 = pdb_origin_graph
    graph2 = graph_pdb_work

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

    result_ZD = mean_squared_error(get_zd_pdb_db(pdb_name, resolution), get_zd_pdb_db(pdb_work, resolution))
    data_write = [[pdb_name, pdb_work, get_chains_pdb_db(pdb_name), get_chains_pdb_db(pdb_work),
                   center_point1, center_point2,
                   center_point1_1, center_point2_1,
                   distance_3d_points(center_point1, center_point1_1),
                   distance_3d_points(center_point2, center_point2_1),
                   resolution, score_pdb_work, result_ZD,
                   time_center, time_center_test,
                   time_aligning]]

    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, data_write)
  if len(list_possibles_pdb) == 0:
    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, [[]])
