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


def do_parallel_test_a(path_data, result_cvs_chain, result_cvs_struct, resolution_range=[5.0, 5.0], can_elements=None,
                       ignore_pdbs=[], percentage_data_set=10, file_checkpoint='check_expe_1a.pkl',
                       error_file='error.txt', can_chain_test=3, can_struct_test=3, add_to_ignore_files=False):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      if not os.path.exists(file_checkpoint):
        all_names = pdb_percentage(percentage_data_set, executor)  # 169315
        open_file = open(file_checkpoint, "wb")
        pickle.dump(all_names, open_file)
        open_file.close()
      else:
        open_file = open(file_checkpoint, "rb")
        all_names = pickle.load(open_file)
        open_file.close()

      # print(all_names, flush=True)
      # all_names = ['1aig']
      # all_names = ['1q8l']
      # all_names = ['7jg5']
      print("Before get pdb names")

      path = os.path.abspath(path_data)
      # print(path)

      if not os.path.isdir(path):
        os.mkdir(path)

      complete_pdb = remove_get_dirs(path_data, can_csv=2, add_to_ignore_files=add_to_ignore_files)
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
        parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_a_aux, path, pdb_name, result_cvs_chain,
                                                        result_cvs_struct, resolution, can_chain_test, can_struct_test),
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


def do_parallel_test_a_aux(path, pdb_name, result_cvs_chain, result_cvs_struct, resolution,
                           can_chain_test, can_struct_test):
  local_path = path + "/" + pdb_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(local_path, pdb_name))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(local_path, pdb_name))

  # print(chains)
  # combinations = combinations_12n(len(chains))[1:]

  start_time = time.time()
  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, pdb_name), path, chains, len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, pdb_name))
  time_eman = time.time() - start_time

  chains_to_segment = {}
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_name, chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    chains_to_segment[chain] = segment
    con_id_segment += 1
  time_segment = time.time() - start_time

  headers_chain = ['Pdb', 'Pdb work', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                   'Point Original syn dis', 'Point Test syn dis',
                   'Resolution',
                   'Score pdb work', 'Changed chain', 'Score chain changed',
                   'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  headers_struct = ['Pdb', 'Pdb work', 'Chains', 'Work Chains', 'Point Original', 'Point Test', 'Point Original syn',
                    'Point Test syn',
                    'Point Original syn dis', 'Point Test syn dis',
                    'Resolution',
                    'Score pdb work',
                    'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  do_test_a_chain(pdb_name, headers_chain, result_cvs_chain, chains_to_segment, resolution, local_path, time_eman,
                  time_segment, can_chain_test)

  do_test_a_struct(pdb_name, headers_struct, result_cvs_struct, chains_to_segment, resolution, local_path, time_eman,
                   time_segment, can_struct_test)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_test_a_struct(pdb_name, headers_csv, result_cvs_file, chains_to_segment, resolution, path_write, time_eman,
                     time_segment, can_chain_test):
  original_segments = []
  for i in chains_to_segment.keys():
    original_segments.append(chains_to_segment[i])

  list_possibles_pdb = get_similar_pdb_struct(pdb_name)

  for i in list_possibles_pdb:
    if i[0] == pdb_name:
      list_possibles_pdb.remove(i)
      break

  random.shuffle(list_possibles_pdb)

  for _i in range(min(can_chain_test, len(list_possibles_pdb))):

    # Change for real method
    pdb_work, score_pdb_work, test_segments, work_chains = get_segments_struct_test(list_possibles_pdb[0], path_write,
                                                                                    resolution)
    list_possibles_pdb.remove(list_possibles_pdb[0])

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

    data_write = [[pdb_name, pdb_work, list(chains_to_segment.keys()), work_chains,
                   center_point1, center_point2,
                   center_point1_1, center_point2_1,
                   distance_3d_points(center_point1, center_point1_1),
                   distance_3d_points(center_point2, center_point2_1),
                   resolution, score_pdb_work,
                   time_segment,
                   time_center, time_center_test,
                   time_graph,
                   time_aligning, time_eman]]

    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, data_write)
  if len(list_possibles_pdb) == 0:
    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, [[]])

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


def do_test_a_chain(pdb_name, headers_csv, result_cvs_file, chains_to_segment, resolution, path_write, time_eman,
                    time_segment, can_chain_test):
  original_segments = []
  list_possibles_pdb_chain = []
  for i in chains_to_segment.keys():
    original_segments.append(chains_to_segment[i])

    # Chains to change
    list_possibles_pdb_chain += get_segments_chain_test(pdb_name, i)

  random.shuffle(list_possibles_pdb_chain)

  for _i in range(min(can_chain_test, len(list_possibles_pdb_chain))):
    pdb_work, score_chain_work, chain_work = list_possibles_pdb_chain.pop()
    id_work = chains_to_segment[chain_work].id_segment

    # Change for real method
    zernike_zd_descriptors_chain_pdb, score_chain_changed = get_zd_descriptors_chain_pdb(pdb_work,
                                                                    chains_to_segment[chain_work].zd_descriptors,
                                                                    path_write, resolution)

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

    data_write = [[pdb_name, pdb_work, list(chains_to_segment.keys()),
                   center_point1, center_point2,
                   center_point1_1, center_point2_1,
                   distance_3d_points(center_point1, center_point1_1),
                   distance_3d_points(center_point2, center_point2_1),
                   resolution, score_chain_work, chain_work, score_chain_changed,
                   time_segment,
                   time_center, time_center_test,
                   time_graph,
                   time_aligning, time_eman]]

    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, data_write)
  if len(list_possibles_pdb_chain) == 0:
    write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, [[]])

def get_segments_chain_test(pdb_name, chain_for_test):
  list_possibles_pdb = get_similar_pdb_chain(pdb_name, chain_for_test)

  for i in list_possibles_pdb:
    if i[0] == pdb_name:
      list_possibles_pdb.remove(i)

  return [i + [chain_for_test] for i in list_possibles_pdb]


def get_zd_descriptors_chain_pdb(pdb_work, zd_descriptors_compare, path_work, resolution):
  download_pdb(pdb_work, '{0}/{1}.pdb'.format(path_work, pdb_work))
  # Maps creation
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_work, pdb_work))

  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(path_work, pdb_work), path_work, chains,
                    len(chains))
  os.remove('{0}/{1}.pdb'.format(path_work, pdb_work))

  best_distance = float('inf')
  result = []

  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one(
      '{0}/{1}_{2}.mrc'.format(path_work + '/' + pdb_work, pdb_work, chain))

    mse = mean_squared_error(zd_descriptors_compare, segments_graph_simulate[0].zd_descriptors)
    if mse < best_distance:
      best_distance = mse
      result = segments_graph_simulate[0].zd_descriptors

  shutil.rmtree(path_work + '/' + pdb_work)

  return result, best_distance
