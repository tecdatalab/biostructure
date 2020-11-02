import os
import pickle

from general_utils.download_utils import get_all_pdb_name, download_pdb, get_all_emd_name, download_emd
from general_utils.list_utils import combinations_12n, get_element_list
from general_utils.math_utils import distance_3d_points
from pdb_to_mrc.miscellaneous import get_chains
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one, get_mrc_segments, get_mrc_synthetic_segments_pdb
from process_mrc.miscellaneous import get_center_point
from writers.csv_writer import write_in_file
import random
import progressbar
import time
import shutil
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

def remove_get_dirs(path):
  result = []
  complete_path = os.path.abspath(path)
  list_dirs = os.listdir(complete_path)

  for dir_name in list_dirs:
    check_path = '{0}/{1}'.format(complete_path, dir_name)

    files_dir = os.listdir(check_path)

    if len(files_dir) == 1 and files_dir[0].split('.')[1]=='csv':
      result.append(dir_name)
    else:
      shutil.rmtree(check_path)
  return result

def do_parallel_test_a(path_data, result_cvs_file, resolution_range=[5.0, 5.0], can_elements=None, remove_files=True,
                       start=None, ignore_pdbs=[]):
  all_names = get_all_pdb_name()
  complete_pdb = remove_get_dirs(path_data)
  ignore_pdbs= ignore_pdbs + complete_pdb
  # all_names = ['101m']
  print("Before get pdb names")

  path = '{0}'.format(os.path.abspath(path_data))
  # print(path)

  if not os.path.isdir(path):
    os.mkdir(path)

  if can_elements is None:
    can_elements = len(all_names)

  bar = progressbar.ProgressBar(maxval=can_elements)
  bar.start()
  con = 0
  con2 = 0
  flag = False
  if start == None:
    flag = True

  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:
      futures = []

      for pdb_name in all_names[:can_elements]:
        if flag == False:
          if pdb_name == start:
            flag = True
          else:
            con += 1
            con2 += 1
            continue
        if pdb_name in ignore_pdbs:
          con += 1
          con2 += 1
          continue

        resolution = random.uniform(resolution_range[0], resolution_range[1])
        # resolution = 3.8680

        con2 += 1
        #print(pdb_name, con2/can_elements)
        futures.append([pdb_name, executor.submit(do_parallel_test_a_aux, path, pdb_name, result_cvs_file,
                                                  remove_files, resolution), resolution])

      for f in futures:
        try:
          d = f[1].result()
          print(f[0], d)
        except Exception as e:
          with open("error_log.txt", "a+") as myfile:
            myfile.write(f[0])
            myfile.write("\n")
            myfile.write(str(f[2]))
            myfile.write("\n")
            myfile.write(str(e))
            myfile.write("\n\n\n\n")
        con += 1
        bar.update(con)


def do_parallel_test_a_aux(path, pdb_name, result_cvs_file, remove_files, resolution):
  local_path = path + "/" + pdb_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(local_path, pdb_name))
  # Maps creation
  chains = get_chains('{0}/{1}.pdb'.format(local_path, pdb_name))
  # print(chains)
  combinations = combinations_12n(len(chains))[1:]

  pdb_to_mrc_chains(True, False, resolution, '{0}/{1}.pdb'.format(local_path, pdb_name), path)
  import sys
  sys.stdout.flush()

  with open('{0}/all_pdb.blist'.format(local_path), 'wb') as fp:
    list_write = ['{0}_{1}.mrc'.format(pdb_name, i) for i in chains]
    pickle.dump(list_write, fp)

  experiments = []

  for test_combination in combinations:
    chains_use = [chains[x] for x in test_combination]
    pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(local_path, pdb_name), path, chains_use)
    experiments.append('{0}_{1}.mrc'.format(pdb_name, "".join(chains_use)))
    # [0] target , [0] original

  with open('{0}/experiments_pdb.blist'.format(local_path), 'wb') as fp:
    pickle.dump(experiments, fp)
  os.remove('{0}/{1}.pdb'.format(local_path, pdb_name))

  headers_csv = ['pdb', 'pdb test', 'Point Original', 'Point Test', 'Point Original sim', 'Point Test sim',
                 'Point Original syn', 'Point Test syn', 'Point Original sim dis', 'Point Test sim dis',
                 'Point Original syn dis', 'Point Test syn dis', 'resolution']

  do_test_a_aux(local_path, pdb_name, headers_csv, result_cvs_file, remove_files, resolution)


def do_test_a_aux(path_data, pdb_name, headers_csv, result_cvs_file, remove_files, resolution):
  segments_graph_simulate, _ = \
    get_mrc_segments('{0}/{1}.mrc'.format(path_data, pdb_name), 3, 1)

  segments_graph_synthetic = None

  with open('{0}/all_pdb.blist'.format(path_data), 'rb') as fp:
    list_all_pdb = pickle.load(fp)
    segments_graph_synthetic, _ = \
      get_mrc_synthetic_segments_pdb('{0}/{1}.mrc'.format(path_data, pdb_name),
                                     "{0}/".format(path_data), list_all_pdb)

  experiments_file = open('{0}/experiments_pdb.blist'.format(path_data), 'rb')
  experiments_list = pickle.load(experiments_file)
  experiments_file.close()

  for experiment in experiments_list:
    # Generate target points
    start_time = time.time()
    segments_graph_complete_target, _ = \
      get_mrc_one('{0}/{1}'.format(path_data, experiment))
    print("--- %s Tiempo de segmentacion ---" % (time.time() - start_time))

    start_time = time.time()
    graph1_match_index = get_element_list(0, [[1, 1]])
    graph2_match_index = get_element_list(1, [[1, 1]])

    center_point1 = get_center_point(graph1_match_index, segments_graph_complete_target, 0)
    center_point2 = get_center_point(graph2_match_index, segments_graph_complete_target, 0)
    print("--- %s Tiempo de puntos centrales ---" % (time.time() - start_time))

    # print("Point Original: ", center_point1, "Point Test: ", center_point2)

    # Generate data simulate
    start_time = time.time()
    segments_graph_simulate_target, _ = \
      get_mrc_segments('{0}/{1}'.format(path_data, experiment), 3, 1)
    print("--- %s Tiempo de segmentacion ---" % (time.time() - start_time))

    # Generate test simulate
    start_time = time.time()
    graph1 = generate_graph(segments_graph_simulate, 50, 0, 6, 1)
    graph2 = generate_graph(segments_graph_simulate_target, 50, 0, 6, 1)
    print("--- %s Tiempo de generacion de grafos ---" % (time.time() - start_time))
    start_time = time.time()
    result = graph_aligning(graph1, graph2, 1, False)
    print("--- %s Tiempo de alinacion ---" % (time.time() - start_time))

    if result != []:
      graph1_match_index = get_element_list(0, result)
      graph2_match_index = get_element_list(1, result)

      start_time = time.time()
      center_point1_1 = get_center_point(graph1_match_index, segments_graph_simulate, 0)
      center_point2_1 = get_center_point(graph2_match_index, segments_graph_simulate_target, 0)
      print("--- %s Tiempo de puntos centrales ---" % (time.time() - start_time))
    else:
      center_point1_1 = [-1, -1, -1]
      center_point2_1 = [-1, -1, -1]

    # print("Point Original sim: ", center_point1_1, "Point Test sim: ", center_point2_1)

    # Generate test synthetic
    start_time = time.time()
    graph1 = generate_graph(segments_graph_synthetic, 50, 0, 6, 1)
    graph2 = generate_graph(segments_graph_simulate_target, 50, 0, 6, 1)
    print("--- %s Tiempo de generacion de grafos ---" % (time.time() - start_time))
    start_time = time.time()
    result = graph_aligning(graph1, graph2, 1, False)
    print("--- %s Tiempo de alineacin ---" % (time.time() - start_time))

    if result != []:

      graph1_match_index = get_element_list(0, result)
      graph2_match_index = get_element_list(1, result)

      start_time = time.time()
      center_point1_2 = get_center_point(graph1_match_index, segments_graph_synthetic, 0)
      center_point2_2 = get_center_point(graph2_match_index, segments_graph_simulate_target, 0)
      print("--- %s Tiempo de puntos centrales ---" % (time.time() - start_time))
    else:
      center_point1_2 = [-1, -1, -1]
      center_point2_2 = [-1, -1, -1]

    # print("Point Original syn: ", center_point1, "Point Test syn: ", center_point2)

    data_write = [[pdb_name, experiment, center_point1, center_point2, center_point1_1, center_point2_1,
                   center_point1_2, center_point2_2,
                   distance_3d_points(center_point1, center_point1_1),
                   distance_3d_points(center_point2, center_point2_1),
                   distance_3d_points(center_point1, center_point1_2),
                   distance_3d_points(center_point2, center_point2_2),
                   resolution]]
    # print(data_write)

    write_in_file('{0}/{1}'.format(path_data, result_cvs_file), headers_csv, data_write)

  if remove_files:
    dirs = os.listdir(path_data)

    for directory in dirs:
      if directory.split('.')[1] != 'csv':
        path_remove = '{0}/{1}'.format(path_data, directory)
        os.remove(path_remove)


def do_parallel_test_b(path_data, result_cvs_file, can_elements=None, remove_files=True):
  all_names = get_all_emd_name()

  path = '{0}'.format(os.path.abspath(path_data))
  # print(path)

  if not os.path.isdir(path):
    os.mkdir(path)

  if can_elements is None:
    can_elements = len(all_names)

  for emd_name in all_names[:can_elements]:
    do_parallel_test_b_aux(path, emd_name, result_cvs_file, remove_files)


def do_parallel_test_b_aux(path, emd_name, result_cvs_file, remove_files):
  local_path = path + "/" + emd_name
  if not os.path.exists(local_path):
    os.makedirs(local_path)

  download_emd(emd_name, '{0}/{1}.map'.format(local_path, emd_name))

  experiments = []
  experiments.append('{0}.map'.format(emd_name))
  # [0] target , [0] original

  with open('{0}/experiments_pdb.blist'.format(local_path), 'wb') as fp:
    pickle.dump(experiments, fp)

  headers_csv = ['emd', 'emd test', 'Point Original', 'Point Test', 'Point Original sim', 'Point Test sim',
                 'Point Original sim dis', 'Point Test sim dis']

  path = '{0}'.format(os.path.abspath(path))
  dirs = os.listdir(path)

  do_test_b_aux(local_path, emd_name, headers_csv, result_cvs_file, remove_files)


def do_test_b_aux(path_data, emd_name, headers_csv, result_cvs_file, remove_files):
  segments_graph_simulate, _ = \
    get_mrc_segments('{0}/{1}.map'.format(path_data, emd_name), 3, 1)
  print('{0}/{1}.map'.format(path_data, emd_name))
  print(len(segments_graph_simulate))

  experiments_file = open('{0}/experiments_pdb.blist'.format(path_data), 'rb')
  experiments_list = pickle.load(experiments_file)
  experiments_file.close()

  for experiment in experiments_list:
    # Generate target points
    print("Antes del completo")
    segments_graph_complete_target, _ = \
      get_mrc_one('{0}/{1}'.format(path_data, experiment))

    graph1_match_index = get_element_list(0, [[1, 1]])
    graph2_match_index = get_element_list(1, [[1, 1]])

    center_point1 = get_center_point(graph1_match_index, segments_graph_complete_target, 0)
    center_point2 = get_center_point(graph2_match_index, segments_graph_complete_target, 0)

    # print("Point Original: ", center_point1, "Point Test: ", center_point2)

    # Generate data simulate
    print("Antes del segundo simulado")
    segments_graph_simulate_target, _ = \
      get_mrc_segments('{0}/{1}'.format(path_data, experiment), 3, 1)

    # Generate test simulate
    graph1 = generate_graph(segments_graph_simulate, 50, 0, 6, 1)
    graph2 = generate_graph(segments_graph_simulate_target, 50, 0, 6, 1)
    result = graph_aligning(graph1, graph2, 1, False)
    if result == []:

      graph1_match_index = get_element_list(0, result)
      graph2_match_index = get_element_list(1, result)

      center_point1_1 = get_center_point(graph1_match_index, segments_graph_simulate, 0)
      center_point2_1 = get_center_point(graph2_match_index, segments_graph_simulate_target, 0)

    else:
      center_point1_1 = [-1, -1, -1]
      center_point2_1 = [-1, -1, -1]

    # print("Point Original sim: ", center_point1_1, "Point Test sim: ", center_point2_1)

    data_write = [[emd_name, experiment, center_point1, center_point2, center_point1_1, center_point2_1,
                   distance_3d_points(center_point1, center_point1_1),
                   distance_3d_points(center_point2, center_point2_1)]]
    # print(data_write)

    write_in_file('{0}/{1}'.format('{0}'.format(path_data), result_cvs_file), headers_csv, data_write)

  if remove_files:
    dirs = os.listdir(path_data)

    for directory in dirs:
      if directory.split('.')[1] != 'csv':
        path_remove = '{0}/{1}'.format(path_data, directory)
        os.remove(path_remove)
