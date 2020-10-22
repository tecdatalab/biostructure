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


def generate_test_data_a(path_data, resolution, can_elements=None):
    all_names = get_all_pdb_name()

    path = '{0}'.format(os.path.abspath(path_data))
    # print(path)
    if not os.path.isdir(path):
        os.mkdir(path)

    if can_elements is None:
        can_elements = len(all_names)

    for pdb_name in all_names[:can_elements]:
        download_pdb(pdb_name, '{0}/{1}.pdb'.format(path, pdb_name))
        # Maps creation
        chains = get_chains('{0}/{1}.pdb'.format(path, pdb_name))
        # print(chains)
        combinations = combinations_12n(len(chains))[1:]

        pdb_to_mrc_chains(True, False, resolution, '{0}/{1}.pdb'.format(path, pdb_name), path)

        with open('{0}/{1}/all_pdb.blist'.format(path, pdb_name), 'wb') as fp:
            list_write = ['{0}_{1}.mrc'.format(pdb_name, i) for i in chains]
            pickle.dump(list_write, fp)

        experiments = []

        for test_combination in combinations:
            chains_use = [chains[x] for x in test_combination]
            pdb_to_mrc_chains(False, False, resolution, '{0}/{1}.pdb'.format(path, pdb_name), path, chains_use)
            experiments.append('{0}_{1}.mrc'.format(pdb_name, "".join(chains_use)))
            # [0] target , [0] original

        with open('{0}/{1}/experiments_pdb.blist'.format(path, pdb_name), 'wb') as fp:
            pickle.dump(experiments, fp)
        os.remove('{0}/{1}.pdb'.format(path, pdb_name))


def do_test_a(path_data, result_cvs_file):
    headers_csv = ['pdb', 'pdb test', 'Point Original', 'Point Test', 'Point Original sim', 'Point Test sim',
                   'Point Original syn', 'Point Test syn', 'Point Original sim dis', 'Point Test sim dis',
                   'Point Original syn dis', 'Point Test syn dis']

    path = '{0}'.format(os.path.abspath(path_data))
    dirs = os.listdir(path)

    for directory in dirs:
        do_test_a_aux(path, directory, headers_csv, result_cvs_file)


def do_test_a_aux(path_data, pdb_name, headers_csv, result_cvs_file):
    segments_graph_simulate, _ = \
        get_mrc_segments('{0}/{1}/{1}.mrc'.format(path_data, pdb_name), 3, 1)

    segments_graph_synthetic = None

    with open('{0}/{1}/all_pdb.blist'.format(path_data, pdb_name), 'rb') as fp:
        list_all_pdb = pickle.load(fp)
        segments_graph_synthetic, _ = \
            get_mrc_synthetic_segments_pdb('{0}/{1}/{1}.mrc'.format(path_data, pdb_name),
                                           "{0}/{1}/".format(path_data, pdb_name), list_all_pdb)

    experiments_file = open('{0}/{1}/experiments_pdb.blist'.format(path_data, pdb_name), 'rb')
    experiments_list = pickle.load(experiments_file)
    experiments_file.close()

    for experiment in experiments_list:
        # Generate target points
        segments_graph_complete_target, _ = \
            get_mrc_one('{0}/{1}/{2}'.format(path_data, pdb_name, experiment))

        graph1_match_index = get_element_list(0, [[1, 1]])
        graph2_match_index = get_element_list(1, [[1, 1]])

        center_point1 = get_center_point(graph1_match_index, segments_graph_complete_target, 0)
        center_point2 = get_center_point(graph2_match_index, segments_graph_complete_target, 0)

        # print("Point Original: ", center_point1, "Point Test: ", center_point2)

        # Generate data simulate
        segments_graph_simulate_target, _ = \
            get_mrc_segments('{0}/{1}/{2}'.format(path_data, pdb_name, experiment), 3, 1)

        # Generate test simulate
        graph1 = generate_graph(segments_graph_simulate, 50, 0, 6, 1)
        graph2 = generate_graph(segments_graph_simulate_target, 50, 0, 6, 1)
        result = graph_aligning(graph1, graph2, 1, False)

        if result != []:
          graph1_match_index = get_element_list(0, result)
          graph2_match_index = get_element_list(1, result)

          center_point1_1 = get_center_point(graph1_match_index, segments_graph_simulate, 0)
          center_point2_1 = get_center_point(graph2_match_index, segments_graph_simulate_target, 0)
        else:
          center_point1_1 = [-1, -1, -1]
          center_point2_1 = [-1, -1, -1]

        # print("Point Original sim: ", center_point1_1, "Point Test sim: ", center_point2_1)

        # Generate test synthetic
        graph1 = generate_graph(segments_graph_synthetic, 50, 0, 6, 1)
        graph2 = generate_graph(segments_graph_simulate_target, 50, 0, 6, 1)
        result = graph_aligning(graph1, graph2, 1, False)

        if result != []:

          graph1_match_index = get_element_list(0, result)
          graph2_match_index = get_element_list(1, result)

          center_point1_2 = get_center_point(graph1_match_index, segments_graph_synthetic, 0)
          center_point2_2 = get_center_point(graph2_match_index, segments_graph_simulate_target, 0)
        else:
          center_point1_1 = [-1, -1, -1]
          center_point2_1 = [-1, -1, -1]

        # print("Point Original syn: ", center_point1, "Point Test syn: ", center_point2)

        data_write = [[pdb_name, experiment, center_point1, center_point2, center_point1_1, center_point2_1,
                       center_point1_2, center_point2_2,
                       distance_3d_points(center_point1, center_point1_1),
                       distance_3d_points(center_point2, center_point2_1),
                       distance_3d_points(center_point1, center_point1_2),
                       distance_3d_points(center_point2, center_point2_2)]]
        # print(data_write)

        write_in_file('{0}/{1}'.format('{0}/{1}'.format(path_data, pdb_name), result_cvs_file), headers_csv, data_write)


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
        segments_graph_complete_target, _ = \
            get_mrc_one('{0}/{1}'.format(path_data, experiment))

        graph1_match_index = get_element_list(0, [[1, 1]])
        graph2_match_index = get_element_list(1, [[1, 1]])

        center_point1 = get_center_point(graph1_match_index, segments_graph_complete_target, 0)
        center_point2 = get_center_point(graph2_match_index, segments_graph_complete_target, 0)

        # print("Point Original: ", center_point1, "Point Test: ", center_point2)

        # Generate data simulate
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
        if directory.split('.')[1] != '.cvs':
          os.remove(directory)
