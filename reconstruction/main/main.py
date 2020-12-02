import os
import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from animation.experiment1_animation import create_ani_expe_1a
from csv_modules.csv_combine import combine_files
from fit.fit_map_chimeraX import fit_map_in_map
from reconstruction.semi_exact_cover import get_semi_exact_s
from csv_modules.csv_writer import write_in_file
from general_utils.list_utils import get_element_list, generate_binary_matrix
from general_utils.math_utils import chance_base_point, get_vector_move_1to2
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from pdb_to_mrc.miscellaneous import get_chains, move_pdb_center
from general_utils.download_utils import download_pdb, download_emd, get_all_pdb_name
from process_graph.graph_algorithm import graph_aligning
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_segments, \
  get_mrc_synthetic_segments_pdb, get_mrc_one
from process_mrc.miscellaneous import get_center_point, \
  get_cube_len_angstrom, get_mass_angstrom, get_mrc_level
from globals.global_values import maps_with_pdb_origin, maps_with_pdb_origin_problems
from metric.metrics_mrc import get_geometric_overlap_p, get_cross_correlation


# segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
#
# segments_graph1 = segments
# segments_graph2 = segments
#
# print(segments[0].mask.shape)
# ##segments = get_mrc_synthetic_segments_pdb("../pdb_to_mrc/exit_pdb/175d", 7)
# graph1 = generate_graph(segments_graph1, 50, 0, 6, 1) #Preguntar con los no conectados y sub grafos
# graph2 = generate_graph(segments_graph2, 50, 0, 6, 1) #Preguntar con los no conectados y sub grafos
#
# result = graph_aligning(graph1, graph2, 1, False)
# print(result)
# #draw_graph_similarity(graph, graph, result)
# #draw_graph_similarity_same_image(graph, graph, result)
# figure_math1,figure_math2 = get_similarity_complete_cloud(segments_graph1,segments_graph2,result)
#
#
# from functools import partial
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from cloud_point import AffineRegistration
# import numpy as np
# import matplotlib
#
# def visualize(iteration, error, X, Y, ax):
#     plt.cla()
#     ax.scatter(X[:, 0],  X[:, 1], X[:, 2], color='red', label='Target')
#     ax.scatter(Y[:, 0],  Y[:, 1], Y[:, 2], color='blue', label='Source')
#     ax.text2D(0.87, 0.92, 'Iteration: {:d}\nQ: {:06.4f}'.format(
#         iteration, error), horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize='x-large')
#     ax.legend(loc='upper left', fontsize='x-large')
#     plt.draw()
#     plt.pause(0.001)
#
#
# def main():
#     matplotlib.use('TKAgg')
#     X = figure_math1
#
#     # synthetic data, equaivalent to X + 1
#     Y = figure_math2
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     callback = partial(visualize, ax=ax)
#
#     print(X.shape)
#     print(Y.shape)
#
#     reg = AffineRegistration(**{'X': X, 'Y': Y})
#     reg.register(callback)
#     plt.show()

def fitting_process(verbose, path_map1, path_map2, segments_graph1, segments_graph2, figure1_shape, figure2_shape,
                    n_points_face, filter_value, max_distance, min_points, path_exit_folder, attempts):
  if verbose:
    print("Fit process start\n\n")

  # Graphs generation
  try:
    graph1 = generate_graph(segments_graph1, n_points_face, filter_value, max_distance, min_points)
    graph2 = generate_graph(segments_graph2, n_points_face, filter_value, max_distance, min_points)
  except Exception as e:
    raise Exception('Error to generate graph, the error is : {0}'.format(str(e)))

  # Match graphs
  try:
    result = graph_aligning(graph1, graph2, 1, False)
  except Exception as e:
    raise Exception('Error to graph aligning, the error is : {0}'.format(str(e)))

  if verbose:
    print("Result match graph", result)

  # Get match segments index
  try:
    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)
  except Exception as e:
    raise Exception('Error to generate graph indexs, the error is : {0}'.format(str(e)))

  # Centers points
  try:
    center_point1 = get_center_point(graph1_match_index, segments_graph1, 0)
    center_point2 = get_center_point(graph2_match_index, segments_graph2, 0)
  except Exception as e:
    raise Exception('Error to calculate centers points, the error is : {0}'.format(str(e)))

  return fitting_process_aux(attempts, center_point1, center_point2, figure1_shape, figure2_shape, path_exit_folder,
                             path_map1, path_map2, verbose)


def fitting_process_all_file(verbose, path_map1, path_map2, segments_graph1, segments_graph2, figure1_shape,
                             figure2_shape, path_exit_folder, attempts):
  if verbose:
    print("Fit process start\n\n")

  # Get match segments index
  try:
    graph1_match_index = get_element_list(0, [[1, 1]])
    graph2_match_index = get_element_list(1, [[1, 1]])
  except Exception as e:
    raise Exception('Error to generate graph indexs, the error is : {0}'.format(str(e)))

  # Centers points
  try:
    center_point1 = get_center_point(graph1_match_index, segments_graph1, 0)
    center_point2 = get_center_point(graph2_match_index, segments_graph2, 0)
  except Exception as e:
    raise Exception('Error to calculate centers points, the error is : {0}'.format(str(e)))

  return fitting_process_aux(attempts, center_point1, center_point2, figure1_shape, figure2_shape, path_exit_folder,
                             path_map1, path_map2, verbose)


def fitting_process_aux(attempts, center_point1_o, center_point2_o, figure1_shape, figure2_shape, path_exit_folder,
                        path_map1, path_map2, verbose):
  center_point1 = center_point1_o.copy()
  center_point2 = center_point2_o.copy()

  if verbose:
    print("Center point map1 grid:", center_point1)
    print("Center point map2 grid:", center_point2)

  # Repair error in 'y' dimension on grid of chimeraX
  center_point1[1] = figure1_shape[1] - center_point1[1]
  center_point2[1] = figure2_shape[1] - center_point2[1]

  if verbose:
    print("Center point map1 grid with y ok:", center_point1)
    print("Center point map2 grid with y ok:", center_point2)

  # Get Angstrom shape
  try:
    real_shape_cube1 = get_cube_len_angstrom(path_map1)
    real_shape_cube2 = get_cube_len_angstrom(path_map2)
  except Exception as e:
    raise Exception('Error to calculate real shape cube, the error is : {0}'.format(str(e)))

  if verbose:
    print("Angstrom shape map1:", real_shape_cube1)
    print("Angstrom shape map2:", real_shape_cube2)

  # Transformation by a rule of 3 from the central point of the grid to the Angstroms cube
  try:
    center_point1_a = chance_base_point(center_point1, figure1_shape, real_shape_cube1)
    center_point2_a = chance_base_point(center_point2, figure2_shape, real_shape_cube2)
  except Exception as e:
    raise Exception('Error to transform center points by rule of 3, the error is : {0}'.format(str(e)))

  if verbose:
    print("Center point in Angstrom cube of map1:", center_point1_a)
    print("Center point in Angstrom cube of map2:", center_point2_a)

  # Obtain the move vector from a 3d point to another 3d point
  try:
    move_vector = get_vector_move_1to2(center_point1_a, center_point2_a)
  except Exception as e:
    raise Exception('Error to calculate move vector, the error is : {0}'.format(str(e)))

  if verbose:
    print("Move vector of map1:", move_vector)

  # Fit map in map
  try:
    result = fit_map_in_map(path_map1, path_map2, path_exit_folder, attempts, map0_vector_move=move_vector)
  except Exception as e:
    raise Exception('Error to calculate fit map in map, the error is : {0}'.format(str(e)))

  result.shape_cube1 = figure1_shape
  result.shape_cube2 = figure2_shape
  result.center_point1_o = center_point1_o
  result.center_point2_o = center_point2_o
  result.center_point1 = center_point1
  result.center_point2 = center_point2
  result.real_shape_cube1 = real_shape_cube1
  result.real_shape_cube2 = real_shape_cube2
  result.center_point1_a = center_point1_a
  result.center_point2_a = center_point2_a
  result.move_vector_map1 = move_vector
  result.percentage_of_overlap = result.overlap_mass / (
    min(get_mass_angstrom(path_map1), get_mass_angstrom(path_map2)))

  if verbose:
    print("Data for result fit are:")
    result.print_data()

  return result


def main_1():
  path = './maps_pdb'
  if not os.path.isdir(path):
    os.mkdir(path)

  download_pdb('175d', '{0}/175d.pdb'.format(path))
  download_pdb('6m03', '{0}/6m03.pdb'.format(path))

  # Maps creation
  chains = get_chains('{0}/175d.pdb'.format(path))
  pdb_to_mrc_chains(True, False, 5.0, '{0}/175d.pdb'.format(path), path, chains, 5)

  chains = get_chains('{0}/6m03.pdb'.format(path))
  pdb_to_mrc_chains(True, False, 5.0, '{0}/6m03.pdb'.format(path), path, chains, 5)

  # Segments generation
  segments_graph1, original_structure1 = get_mrc_segments("{0}/175d/175d.mrc".format(path), 7, 3, 1)
  segments_graph2, original_structure2 = get_mrc_segments("{0}/6m03/6m03.mrc".format(path), 7, 3, 1)

  figure1_shape = original_structure1.mask.shape
  figure2_shape = original_structure2.mask.shape

  result = fitting_process(True, "{0}/175d/175d.mrc".format(path), "{0}/6m03/6m03.mrc".format(path), segments_graph1,
                           segments_graph2, figure1_shape, figure2_shape, 50, 0, 6, 1, './exit_fit', 50)


def main_2():
  path = './maps_pdb'
  if not os.path.isdir(path):
    os.mkdir(path)

  download_pdb('175d', '{0}/175d.pdb'.format(path))

  # Maps creation
  chains = get_chains('{0}/175d.pdb'.format(path))
  pdb_to_mrc_chains(True, False, 5.0, '{0}/175d.pdb'.format(path), path, chains, len(chains))

  # Segments generation
  segments_graph1, original_structure1 = get_mrc_synthetic_segments_pdb("{0}/175d/175d.mrc".format(path),
                                                                        "{0}/175d".format(path), 7)
  figure1_shape = original_structure1.mask.shape

  result = fitting_process(True, "{0}/175d/175d.mrc".format(path), "{0}/175d/175d.mrc".format(path), segments_graph1,
                           segments_graph1, figure1_shape, figure1_shape, 50, 0, 6, 1, './exit_fit', 50)


def main_3():
  con = 0
  len_file = len(maps_with_pdb_origin)
  headers_csv = ['map_name', 'center_p1A', 'center_p2A', 'move_v1', 'move_v2', 'map1_level', 'map2_level',
                 'Num_poins',
                 'Correlation', 'Correlation_about_mean', 'Overlap', 'Steps', 'Shift', 'Angle', 'Matrix_rt', 'Axis',
                 'Axis_point', 'Rotation_angle', 'Shift_along_axis', 'map1_path', 'map2_path']
  initial_flag = False
  for i in maps_with_pdb_origin:
    con += 1
    name = i[4:-1]
    name = name[:-3]
    if name in maps_with_pdb_origin_problems:
      continue
    if name == '0009':
      initial_flag = True
    if not initial_flag:
      continue

    print("Actual execution: ", name)

    # Creation of segments in this case of only 1 element
    try:
      segments_graph1, original_structure1 = \
        get_mrc_one("/mnt/hgfs/Project_files/selected_sim/sim_emd_{0}.mrc".format(name), 7)
      segments_graph2, original_structure2 = \
        get_mrc_one("/mnt/hgfs/Project_files/original_maps/original_maps/emd_{0}.map".format(name), 7)
      figure1_shape = original_structure1.mask.shape
      figure2_shape = original_structure2.mask.shape
    except Exception as e:
      f = open('log.txt', "a+")
      f.write('Error to open file: {0}'.format(name))
      f.write(str(e))
      f.write(str('\n\n'))
      f.close()
      continue

    # Fit map in map
    result = fitting_process_all_file(True, '/mnt/hgfs/Project_files/selected_sim/sim_emd_{0}.mrc'.format(name),
                                      '/mnt/hgfs/Project_files/original_maps/original_maps/emd_{0}.map'.format(
                                        name),
                                      segments_graph1, segments_graph2, figure1_shape, figure2_shape,
                                      './exit_fit_god',
                                      50)

    data_write = [[name, str(result.center_point1_o), str(result.center_point2_o), str(result.move_vector_map1),
                   [0, 0, 0], str('Defa'), str('Defa'), str(result.num_poins), str(result.correlation),
                   str(result.correlation_about_mean), str(result.overlap), str(result.steps), str(result.shift),
                   str(result.angle), str(result.matrix_rt), str(result.axis), str(result.axis_point),
                   str(result.rotation_angle), str(result.shift_along_axis),
                   '/mnt/hgfs/Project_files/selected_sim/sim_emd_{0}.mrc'.format(name),
                   '/mnt/hgfs/Project_files/original_maps/original_maps/emd_{0}.map'.format(name)]]

    write_in_file('./exit_fit_god/results2.csv', headers_csv, data_write)
    result.print_data()
    print("Actual execution: ", con, "  of map: ", name, " progress: ", con / len_file)


def main_4():
  segments_graph1, original_structure1 = get_mrc_one("/home/lcastillo98/Documents/git_projects/sim_emd_9882.mrc", 7)
  segments_graph2, original_structure2 = get_mrc_one("/home/lcastillo98/Documents/git_projects/emd_9882.map", 7)

  figure1_shape = original_structure1.mask.shape
  figure2_shape = original_structure2.mask.shape

  # Fit map in map
  result = fitting_process_all_file(True, "/home/lcastillo98/Documents/git_projects/sim_emd_9882.mrc",
                                    "/home/lcastillo98/Documents/git_projects/emd_9882.map",
                                    segments_graph1, segments_graph2, figure1_shape, figure2_shape, './exit_fit',
                                    50)
  print("Fin")


def main_5():
  download_emd("0009", './exit_fit/0009.map', True)


def main_6():
  import mrcfile

  file1 = mrcfile.open('/home/lcastillo98/Documents/git_projects/sim_emd_9882_fit.mrc')
  file2 = mrcfile.open('/home/lcastillo98/Documents/git_projects/emd_9882.map')

  # file1 = mrcfile.open('/home/lcastillo98/Documents/git_projects/cube50x50x50.map')
  # file2 = mrcfile.open('/home/lcastillo98/Documents/git_projects/cube100x100x100.map')

  percentage_overlap = get_geometric_overlap_p(file1.data, file2.data)
  print("Percentage of overlap:", percentage_overlap)

  cross_correlation = get_cross_correlation(file1.data, file2.data)
  print("Cross_correlation:", cross_correlation)


def main_7():
  initial_matrix = [[1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 0, 0, 1, 0, 0],
                    [0, 0, 1, 1, 0, 1, 1],
                    [0, 0, 0, 1, 0, 0, 1],
                    [1, 0, 1, 0, 0, 1, 0],
                    [1, 1, 1, 1, 1, 1, 1]
                    ]
  top = 12
  binary_matrix = generate_binary_matrix(initial_matrix)
  combinations = get_semi_exact_s(binary_matrix, top, 2)
  print("Combinations: ", combinations)


def main_8():
  path = './maps_pdb'
  if not os.path.isdir(path):
    os.mkdir(path)

  download_pdb('175d', '{0}/175d.pdb'.format(path))

  # Maps creation
  chains = get_chains('{0}/175d.pdb'.format(path))
  pdb_to_mrc_chains(True, False, 5.0, '{0}/175d.pdb'.format(path), path, chains, len(chains))

  # Segments generation
  segments_graph1, original_structure1 = get_mrc_synthetic_segments_pdb("{0}/175d/175d.mrc".format(path),
                                                                        "{0}/175d".format(path), 7)
  initial_matrix = []
  for i in segments_graph1:
    flat_data = i.mask.ravel()
    flat_data[flat_data > 0] = 1
    flat_data = flat_data.astype(int)
    flat_data = flat_data.tolist()

    initial_matrix.append(flat_data)
  print("Can elements ", len(initial_matrix))
  print(initial_matrix)

  top = 10
  binary_matrix = generate_binary_matrix(initial_matrix)
  combinations = get_semi_exact_s(binary_matrix, len(initial_matrix[0]), top, 80)
  print("Combinations: ", combinations)


def main_9():
  import itertools
  import time
  import matplotlib.pyplot as plt
  times = []
  cans = []

  for i in range(2, 11):
    initial_matrix = list(itertools.product([0, 1], repeat=i))
    top = 10000
    binary_matrix = generate_binary_matrix(initial_matrix)
    can_elements = len(initial_matrix)
    start_time = time.time()
    combinations = get_semi_exact_s(binary_matrix, top, 20)
    total_seconds = (time.time() - start_time)
    print("--- {0} seconds --- can elements: {1}".format(total_seconds, can_elements))
    times.append(total_seconds)
    cans.append(can_elements)

  plt.plot(cans, times)
  plt.ylabel('Time in seconds')
  plt.xlabel('Can elements')
  plt.show()


def main_10():
  # path = '/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/data_experiment_1_a'
  #
  # create_ani_expe_1a(path,
  #                    '1bgy',
  #                    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
  #                     'U', 'V', 'W'],
  #                    ['A', 'B', 'C', 'D', 'E'],
  #                    [[1, 1], [1, 2], [2, 3]],
  #                    0.1,
  #                    [98, 98, 98],
  #                    [150, 150, 150],
  #                    [150, 150, 150],
  #                    [150, 150, 150],
  #                    5.0,
  #                    True)

  from pymol import cmd

  from pymol import cmd

  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy.mrc", finish=0)
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_A.mrc", finish=1)
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_B.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_C.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_D.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_E.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_F.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_G.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_H.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_I.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_J.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_K.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_M.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_N.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_O.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_P.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_Q.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_R.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_S.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_T.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_U.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_V.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_W.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_father_pc.mrc")
  cmd.load("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1bgy_father_pt.mrc")

  cmd.volume("1bgy_volume", "1bgy")
  cmd.show("volume", "1bgy_volume")

  cmd.isosurface("1bgy_A_surface", "1bgy_A")
  cmd.hide("surface", "1bgy_A_surface")
  cmd.isosurface("1bgy_B_surface", "1bgy_B")
  cmd.hide("surface", "1bgy_B_surface")
  cmd.isosurface("1bgy_C_surface", "1bgy_C")
  cmd.hide("surface", "1bgy_C_surface")
  cmd.isosurface("1bgy_D_surface", "1bgy_D")
  cmd.hide("surface", "1bgy_D_surface")
  cmd.isosurface("1bgy_E_surface", "1bgy_E")
  cmd.hide("surface", "1bgy_E_surface")
  cmd.isosurface("1bgy_F_surface", "1bgy_F")
  cmd.hide("surface", "1bgy_F_surface")
  cmd.isosurface("1bgy_G_surface", "1bgy_G")
  cmd.hide("surface", "1bgy_G_surface")
  cmd.isosurface("1bgy_H_surface", "1bgy_H")
  cmd.hide("surface", "1bgy_H_surface")
  cmd.isosurface("1bgy_I_surface", "1bgy_I")
  cmd.hide("surface", "1bgy_I_surface")
  cmd.isosurface("1bgy_J_surface", "1bgy_J")
  cmd.hide("surface", "1bgy_J_surface")
  cmd.isosurface("1bgy_K_surface", "1bgy_K")
  cmd.hide("surface", "1bgy_K_surface")
  cmd.isosurface("1bgy_M_surface", "1bgy_M")
  cmd.hide("surface", "1bgy_M_surface")
  cmd.isosurface("1bgy_N_surface", "1bgy_N")
  cmd.hide("surface", "1bgy_N_surface")
  cmd.isosurface("1bgy_O_surface", "1bgy_O")
  cmd.hide("surface", "1bgy_O_surface")
  cmd.isosurface("1bgy_P_surface", "1bgy_P")
  cmd.hide("surface", "1bgy_P_surface")
  cmd.isosurface("1bgy_Q_surface", "1bgy_Q")
  cmd.hide("surface", "1bgy_Q_surface")
  cmd.isosurface("1bgy_R_surface", "1bgy_R")
  cmd.hide("surface", "1bgy_R_surface")
  cmd.isosurface("1bgy_S_surface", "1bgy_S")
  cmd.hide("surface", "1bgy_S_surface")
  cmd.isosurface("1bgy_T_surface", "1bgy_T")
  cmd.hide("surface", "1bgy_T_surface")
  cmd.isosurface("1bgy_U_surface", "1bgy_U")
  cmd.hide("surface", "1bgy_U_surface")
  cmd.isosurface("1bgy_V_surface", "1bgy_V")
  cmd.hide("surface", "1bgy_V_surface")
  cmd.isosurface("1bgy_W_surface", "1bgy_W")
  cmd.hide("surface", "1bgy_W_surface")
  cmd.isosurface("1bgy_father_pc_surface", "1bgy_father_pc")
  cmd.color("red", "1bgy_father_pc_surface")
  cmd.isosurface("1bgy_father_pt_surface", "1bgy_father_pt")
  cmd.color("white", "1bgy_father_pt_surface")

  cmd.png("/home/lcastillo98/Desktop/data_experiment_1_a/1bgy/1.png")


def experiment_1():
  from experiment.experiment_1 import do_parallel_test_a, do_parallel_test_b

  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/data_experiment_1_a".format(local_path), "result.csv", [3.5, 9.5], 1,
                     ignore_pdbs=['105d', '106d', '108d'])
  print("Finish")

def experiment_3():
  from experiment.experiment_3 import do_parallel_test_a

  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/data_experiment_3_a".format(local_path), "result.csv", [3.5, 9.5], 1)
  print("Finish")


def union_test():
  combine_files('salida.csv',
                '/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/data_experiment_1_a')


if __name__ == '__main__':
  # main_1()
  # main_2()
  # main_3()
  # main_4()
  # main_5()
  # main_6()
  # main_7()
  # main_8()
  # main_9()
  # main_10()
  # experiment_1()
  experiment_3()
  # union_test()
