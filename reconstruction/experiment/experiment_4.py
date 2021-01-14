import os
import random
import shutil
import time

import numpy as np
from memory_profiler import memory_usage
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from experiment.utils_general import get_ignore_pdbs
from general_utils.download_utils import download_pdb
from general_utils.list_utils import generate_binary_matrix
from pdb_to_mrc.miscellaneous import get_chains
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_mrc.generate import get_mrc_one
from reconstruction.DLX import solve, gen_y_dicc, gen_x_dicc
from reconstruction.semi_exact_cover import get_semi_exact_s


def remove_get_dirs(path):
  result = []
  complete_path = os.path.abspath(path)
  list_dirs = os.listdir(complete_path)

  for dir_name in list_dirs:
    check_path = '{0}/{1}'.format(complete_path, dir_name)
    files_dir = os.listdir(check_path)

    if len(files_dir) == 1 and files_dir[0].split('.')[1] == 'csv':
      result.append(dir_name)
    else:
      shutil.rmtree(check_path)

  return result


def do_parallel_test_a(path_data, result_cvs_file, resolution_range=[5.0, 5.0], can_elements=None,
                       ignore_pdbs=[]):
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      # all_names = get_all_pdb_name()  # 169315
      all_names = ['100d']
      # all_names = ['7jsh']
      print("Before get pdb names")

      path = os.path.abspath(path_data)
      # print(path)

      if not os.path.isdir(path):
        os.mkdir(path)

      complete_pdb = remove_get_dirs(path_data)
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
        # parallel_jobs.append([pdb_name, executor.submit(do_parallel_test_a_aux, path, pdb_name, result_cvs_file,
        #                                                  resolution), resolution])
        do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution)
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          with open("error_log.txt", "a+") as myfile:
            myfile.write(f[0])
            myfile.write("\n")
            myfile.write(str(f[2]))
            myfile.write("\n")
            myfile.write(str(type(e).__name__))
            myfile.write("\n")
            myfile.write(str(e))
            myfile.write("\n\n\n\n")


def do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution):
  local_path = path + "/" + pdb_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(local_path, pdb_name))
  # Maps creation
  chains = get_chains('{0}/{1}.pdb'.format(local_path, pdb_name))
  # print(chains)

  start_time = time.time()
  pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, pdb_name), path, chains, len(chains))
  os.remove('{0}/{1}.pdb'.format(local_path, pdb_name))
  time_eman = time.time() - start_time

  all_segments = []
  start_time = time.time()
  con_id_segment = 1
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_name, chain), calculate_Z3D=False)
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    all_segments.append(segment)
    con_id_segment += 1
  time_load = time.time() - start_time

  headers_csv = ['Pdb', 'Chains', 'Resolution',
                 'Time load', 'Time EMAN2', 'Time our method', 'Time DLX method',
                 'Our memory', 'DLX memory',
                 'Num update our', 'Num update DLX']

  do_test_a(pdb_name, headers_csv, result_cvs_file, all_segments, resolution, local_path,
            time_eman, time_load, chains)

  dirs = os.listdir(local_path)

  for directory in dirs:
    if directory.split('.')[1] != 'csv':
      path_remove = '{0}/{1}'.format(local_path, directory)
      os.remove(path_remove)


def do_our_method(initial_matrix, count_update=[]):
  # print('Before')
  binary_matrix = generate_binary_matrix(initial_matrix)
  # print('After')
  combinations = get_semi_exact_s(binary_matrix, 1, 1, count_update=count_update)
  return combinations


def do_DLX_method(initial_matrix, count_update=[]):
  Y = gen_y_dicc(initial_matrix)
  X = gen_x_dicc(Y, initial_matrix)
  result = list(solve(X, Y, count_update=count_update))
  return result


def do_test_a(pdb_name, headers_csv, result_cvs_file, all_segments, resolution, local_path,
              time_eman, time_load, chains):
  aux_cube = np.zeros((all_segments[0].mask.shape[0] * all_segments[0].mask.shape[1] * all_segments[0].mask.shape[2],))
  initial_matrix = []

  for i in range(len(all_segments)):
    initial_matrix.append(all_segments[i].mask.ravel())

  for i in range(len(initial_matrix)):
    points_i = np.where(initial_matrix[i] > 0)[0]
    initial_matrix[i][points_i] = 1
    aux_cube[points_i] = 1
    for j in range(i + 1, len(initial_matrix)):
      points_j = np.where(initial_matrix[j] > 0)[0]
      change_points = np.intersect1d(points_i, points_j)
      initial_matrix[j][change_points] = 0

  free_points = np.where(aux_cube == 0)
  initial_matrix[0][free_points] = 1

  # initial_matrix = [[0, 1, 0],
  #                   [1, 0, 0],
  #                   [0, 0, 1]]

  # Our method
  start_time = time.time()
  combinations = do_our_method(initial_matrix)
  our_combination_time = time.time() - start_time
  # print(our_combination_time)
  max_mem_our = max(memory_usage((do_our_method, (initial_matrix,)), interval=.2))
  # print(max_mem_our)

  # DLX method
  start_time = time.time()
  result = do_DLX_method(initial_matrix)
  dlx_combination_time = time.time() - start_time
  # print(dlx_combination_time)
  max_mem_dlx = max(memory_usage((do_DLX_method, (initial_matrix,)), interval=.2))
  # print(max_mem_dlx)

  count_update_DLX = [0]
  do_DLX_method(initial_matrix, count_update_DLX)
  # print(count_update_DLX)

  count_update_our = [0]
  do_our_method(initial_matrix, count_update_our)
  # print(count_update_our)
  # print(result, result)
  # print(combinations[0][1], combinations[0])
  # print(np.setdiff1d(combinations[0][1], result[0]))

  if not (np.setdiff1d(combinations[0][1], result[0]).tolist() == []):
    raise Exception("The two methods do not produce the same result")

  if not len(combinations[0][1]) == len(initial_matrix):
    raise Exception("Can not get exact cover with our method")

  if not len(result[0]) == len(initial_matrix):
    raise Exception("Can not get exact cover with DLX method")

  data_write = [[pdb_name,
                 chains,
                 resolution,
                 time_load, time_eman,
                 our_combination_time, dlx_combination_time,
                 max_mem_our, max_mem_dlx,
                 count_update_our[0], count_update_DLX[0]]]

  write_in_file('{0}/{1}'.format(local_path, result_cvs_file), headers_csv, data_write)
