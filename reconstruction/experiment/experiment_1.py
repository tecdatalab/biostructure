import os
import random
import time
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import traceback

from csv_modules.csv_writer import write_in_file
from experiment.utils_experiment_1 import gen_keys_experiemnts
from experiment.utils_general import remove_get_dirs
from general_utils.pdb_utils import get_ignore_pdbs
from general_utils.download_utils import get_all_pdb_name, download_pdb
from general_utils.list_utils import get_element_list
from general_utils.math_utils import distance_3d_points
from pdb_to_mrc.miscellaneous import get_chains
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from process_mrc.miscellaneous import get_center_point


def do_parallel_test_a(path_data, result_cvs_file, resolution_range=[5.0, 5.0], can_elements=None,
                      ignore_pdbs=[], range_incompleteness=[10.0, 15.0], can_try_experiments=10,
                       add_to_ignore_files=False, error_file='error.txt'):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      all_names = get_all_pdb_name()  # 169315
      # all_names = ['7jsh']
      print("Before get pdb names")

      path = os.path.abspath(path_data)
      # print(path)

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

        resolution = random.uniform(resolution_range[0], resolution_range[1])
        # resolution = 3.8680

        # print(pdb_name, con2/can_elements)
        parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_a_aux, path, pdb_name, result_cvs_file,
                                                        resolution, range_incompleteness, can_try_experiments),
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


def do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution, range_incompleteness, can_try_experiments):
  local_path = path + "/" + pdb_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(local_path, pdb_name))
  # Maps creation
  chains = get_chains('{0}/{1}.pdb'.format(local_path, pdb_name))

  # print(chains)
  # combinations = combinations_12n(len(chains))[1:]

  start_time = time.time()
  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, pdb_name), path, chains, len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, pdb_name))
  time_eman = time.time() - start_time

  segments_to_chains = {}
  all_segments = []
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_name, chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    segments_to_chains[con_id_segment] = chain
    all_segments.append(segment)
    con_id_segment += 1
  time_segment = time.time() - start_time

  headers_csv = ['Pdb', 'Chains', 'Point Original', 'Point Test', 'Point Original syn', 'Point Test syn',
                 'Point Original syn dis', 'Point Test syn dis',
                 'Resolution',
                 'Match', 'Father Chains', 'Test Chains',
                 'Time segment', 'Time center', 'Time test', 'Time graph', 'Time alignment', 'Time EMAN2']

  completed_experiments = []
  for _ in range(can_try_experiments):
    can_try = 10
    while can_try > 0:
      percentage = random.uniform(min(range_incompleteness), max(range_incompleteness))

      point_test = [int(random.uniform(0, all_segments[0].mask.shape[0])),
                    int(random.uniform(0, all_segments[0].mask.shape[1])),
                    int(random.uniform(0, all_segments[0].mask.shape[2]))]
      # print(point_test)
      key_segment_test = gen_keys_experiemnts(all_segments, 50, 0, [int(random.uniform(0, all_segments[0].mask.shape[0])),
                                                                    int(random.uniform(0, all_segments[0].mask.shape[1])),
                                                                    int(random.uniform(0, all_segments[0].mask.shape[2]))],
                                              percentage)
      # print(key_segment_test)
      if key_segment_test not in completed_experiments:
        # print(key_segment_test)
        completed_experiments.append(key_segment_test)
        test_segments = [i for i in all_segments if i.id_segment in key_segment_test]
        test_chains = [segments_to_chains[i] for i in key_segment_test]

        do_test_a(pdb_name, headers_csv, result_cvs_file, all_segments, test_segments, test_chains, resolution,
                  local_path,
                  time_eman, time_segment, chains, test_chains)
        break
      else:
        can_try -=1

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_test_a(pdb_name, headers_csv, result_cvs_file, all_segments, test_segments, test_chains, resolution, path_write,
              time_eman, time_segment, fchains, tchains):
  false_match_list = [[i.id_segment, i.id_segment] for i in test_segments]
  graph1_match_index = get_element_list(0, false_match_list)
  graph2_match_index = get_element_list(1, false_match_list)

  start_time = time.time()
  center_point1 = get_center_point(graph1_match_index, test_segments, 0)
  center_point2 = get_center_point(graph2_match_index, test_segments, 0)
  time_center = time.time() - start_time

  # Generate test syntetic
  start_time = time.time()
  graph1 = generate_graph(all_segments, 50, 0, 6, 1)
  graph2 = generate_graph(test_segments, 50, 0, 6, 1)
  time_graph = time.time() - start_time

  start_time = time.time()
  result = graph_aligning(graph1, graph2, 2, False)
  time_aligning = time.time() - start_time

  time_center_test = 0
  if result != []:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    start_time = time.time()
    center_point1_1 = get_center_point(graph1_match_index, all_segments, 0)
    center_point2_1 = get_center_point(graph2_match_index, test_segments, 0)
    time_center_test = time.time() - start_time
  else:
    center_point1_1 = [-1, -1, -1]
    center_point2_1 = [-1, -1, -1]

  data_write = [[pdb_name, test_chains,
                 center_point1, center_point2,
                 center_point1_1, center_point2_1,
                 distance_3d_points(center_point1, center_point1_1),
                 distance_3d_points(center_point2, center_point2_1),
                 resolution,
                 result, fchains, tchains,
                 time_segment,
                 time_center, time_center_test,
                 time_graph,
                 time_aligning, time_eman]]

  write_in_file('{0}/{1}'.format(path_write, result_cvs_file), headers_csv, data_write)
