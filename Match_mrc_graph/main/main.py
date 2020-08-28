import shutil

from process_mrc.generate import get_mrc_segments, \
    get_mrc_synthetic_segments_pdb, get_mrc_one
from process_graph.process_graph_utils import generate_graph, \
    draw_graph_similarity, draw_graph_similarity_same_image, \
    get_similarity_complete_cloud
from process_graph.graph_algorithm import graph_aligning
from fit.fit_map import fit_map_in_map
from pdb_mrc.pdb_transform import download_pdb, get_chains, pdb_to_mrc_chains
import os

from process_mrc.miscellaneous import get_center_point, get_vector_move_1to2, \
    get_cube_len, chance_basec, get_mass
from utils_match.csv_writer import write_in_file
from utils_match.list_processing import get_element_list

import numpy as np


# segments = get_mrc_segments("../../maps/1010/EMD-1010.map", 7, 3, 1)
# 
# segments_graph1 = segments
# segments_graph2 = segments
# 
# print(segments[0].mask.shape)
# ##segments = get_mrc_synthetic_segments_pdb("../pdb_mrc/exit_pdb/175d", 7)
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
        real_shape_cube1 = get_cube_len(path_map1)
        real_shape_cube2 = get_cube_len(path_map2)
    except Exception as e:
        raise Exception('Error to calculate real shape cube, the error is : {0}'.format(str(e)))

    if verbose:
        print("Angstrom shape map1:", real_shape_cube1)
        print("Angstrom shape map2:", real_shape_cube2)

    # Transformation by a rule of 3 from the central point of the grid to the Angstroms cube
    try:
        center_point1_a = chance_basec(center_point1, figure1_shape, real_shape_cube1)
        center_point2_a = chance_basec(center_point2, figure2_shape, real_shape_cube2)
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
    result.percentage_of_overlap = result.overlap_mass/(min(get_mass(path_map1), get_mass(path_map2)))

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
    segments_graph1, figure1_shape = get_mrc_segments("{0}/175d/175d.mrc".format(path), 7, 3, 1)
    segments_graph2, figure2_shape = get_mrc_segments("{0}/6m03/6m03.mrc".format(path), 7, 3, 1)

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
    segments_graph1, figure1_shape = get_mrc_synthetic_segments_pdb("{0}/175d/175d.mrc".format(path),
                                                                    "{0}/175d".format(path), 7)

    result = fitting_process(True, "{0}/175d/175d.mrc".format(path), "{0}/175d/175d.mrc".format(path), segments_graph1,
                             segments_graph1, figure1_shape, figure1_shape, 50, 0, 6, 1, './exit_fit', 50)


def main_3():
    files_list = ['emd_0009.map', 'emd_0244.map', 'emd_0293.map', 'emd_0355.map', 'emd_0362.map', 'emd_0366.map',
                  'emd_0367.map', 'emd_0404.map', 'emd_0434.map', 'emd_0493.map', 'emd_0646.map', 'emd_0647.map',
                  'emd_0648.map', 'emd_0695.map', 'emd_0728.map', 'emd_0729.map', 'emd_0777.map', 'emd_0778.map',
                  'emd_0790.map', 'emd_0901.map', 'emd_0979.map', 'emd_10132.map', 'emd_10194.map', 'emd_10373.map',
                  'emd_10464.map', 'emd_10699.map', 'emd_1181.map', 'emd_1261.map', 'emd_1263.map', 'emd_1884.map',
                  'emd_1918.map', 'emd_1983.map', 'emd_1988.map', 'emd_1989.map', 'emd_1990.map', 'emd_2005.map',
                  'emd_20267.map', 'emd_20287.map', 'emd_20493.map', 'emd_20511.map', 'emd_2071.map', 'emd_2072.map',
                  'emd_20721.map', 'emd_20813.map', 'emd_20924.map', 'emd_20931.map', 'emd_21121.map', 'emd_21335.map',
                  'emd_21585.map', 'emd_21602.map', 'emd_2233.map', 'emd_2528.map', 'emd_2594.map', 'emd_2595.map',
                  'emd_2605.map', 'emd_2676.map', 'emd_2706.map', 'emd_2751.map', 'emd_2752.map', 'emd_2786.map',
                  'emd_2813.map', 'emd_2845.map', 'emd_2860.map', 'emd_2917.map', 'emd_2925.map', 'emd_30036.map',
                  'emd_3101.map', 'emd_3133.map', 'emd_3164.map', 'emd_3165.map', 'emd_3166.map', 'emd_3167.map',
                  'emd_3168.map', 'emd_3169.map', 'emd_3170.map', 'emd_3235.map', 'emd_3241.map', 'emd_3305.map',
                  'emd_3329.map', 'emd_3352.map', 'emd_3355.map', 'emd_3356.map', 'emd_3433.map', 'emd_3445.map',
                  'emd_3467.map', 'emd_3468.map', 'emd_3469.map', 'emd_3470.map', 'emd_3471.map', 'emd_3472.map',
                  'emd_3473.map', 'emd_3474.map', 'emd_3475.map', 'emd_3477.map', 'emd_3478.map', 'emd_3479.map',
                  'emd_3480.map', 'emd_3481.map', 'emd_3482.map', 'emd_3483.map', 'emd_3484.map', 'emd_3522.map',
                  'emd_3527.map', 'emd_3529.map', 'emd_3530.map', 'emd_3604.map', 'emd_3634.map', 'emd_3669.map',
                  'emd_3696.map', 'emd_3706.map', 'emd_3778.map', 'emd_3803.map', 'emd_3850.map', 'emd_3894.map',
                  'emd_3896.map', 'emd_3954.map', 'emd_4021.map', 'emd_4025.map', 'emd_4057.map', 'emd_4078.map',
                  'emd_4100.map', 'emd_4154.map', 'emd_4318.map', 'emd_4515.map', 'emd_4671.map', 'emd_4680.map',
                  'emd_4707.map', 'emd_4708.map', 'emd_4710.map', 'emd_4742.map', 'emd_4743.map', 'emd_4744.map',
                  'emd_4764.map', 'emd_4918.map', 'emd_4940.map', 'emd_4975.map', 'emd_5002.map', 'emd_5100.map',
                  'emd_5258.map', 'emd_5335.map', 'emd_5391.map', 'emd_5395.map', 'emd_5423.map', 'emd_5607.map',
                  'emd_5609.map', 'emd_5610.map', 'emd_6102.map', 'emd_6284.map', 'emd_6285.map', 'emd_6286.map',
                  'emd_6369.map', 'emd_6476.map', 'emd_6779.map', 'emd_6781.map', 'emd_6782.map', 'emd_6783.map',
                  'emd_6823.map', 'emd_6826.map', 'emd_6850.map', 'emd_7065.map', 'emd_7090.map', 'emd_7305.map',
                  'emd_7328.map', 'emd_7463.map', 'emd_7529.map', 'emd_7896.map', 'emd_7939.map', 'emd_7940.map',
                  'emd_8132.map', 'emd_8133.map', 'emd_8134.map', 'emd_8148.map', 'emd_8149.map', 'emd_8190.map',
                  'emd_8427.map', 'emd_8513.map', 'emd_8539.map', 'emd_8579.map', 'emd_8886.map', 'emd_8898.map',
                  'emd_8981.map', 'emd_9043.map', 'emd_9134.map', 'emd_9147.map', 'emd_9163.map', 'emd_9272.map',
                  'emd_9302.map', 'emd_9303.map', 'emd_9378.map', 'emd_9387.map', 'emd_9388.map', 'emd_9594.map',
                  'emd_9621.map', 'emd_9714.map', 'emd_9715.map', 'emd_9716.map', 'emd_9757.map', 'emd_9841.map',
                  'emd_9881.map', 'emd_9882.map', 'emd_10944.map', 'emd_10960.map', 'emd_10748.map', 'emd_10754.map',
                  'emd_11173.map']
    files_problems_list = ["3235", "3706", "4680"]
    con = 0
    len_file = len(files_list)
    headers_csv = ['map_name', 'center_p1A', 'center_p2A', 'move_v1', 'move_v2', 'map1_level', 'map2_level',
                   'Num_poins',
                   'Correlation', 'Correlation_about_mean', 'Overlap', 'Steps', 'Shift', 'Angle', 'Matrix_rt', 'Axis',
                   'Axis_point', 'Rotation_angle', 'Shift_along_axis', 'map1_path', 'map2_path']
    initial_flag = False
    for i in files_list:
        con += 1
        name = i[4:-1]
        name = name[:-3]
        if name in files_problems_list:
            continue
        if name == '0009':
            initial_flag = True
        if not initial_flag:
            continue

        print("Actual execution: ", name)

        # Creation of segments in this case of only 1 element
        try:
            segments_graph1, figure1_shape = \
                get_mrc_one("/mnt/hgfs/Project_files/selected_sim/sim_emd_{0}.mrc".format(name), 7)
            segments_graph2, figure2_shape = \
                get_mrc_one("/mnt/hgfs/Project_files/original_maps/original_maps/emd_{0}.map".format(name), 7)
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
    segments_graph1, figure1_shape = get_mrc_one("/home/lcastillo98/Documents/git_projects/sim_emd_9882.mrc", 7)
    segments_graph2, figure2_shape = get_mrc_one("/home/lcastillo98/Documents/git_projects/emd_9882.map", 7)

    # Fit map in map
    result = fitting_process_all_file(True, "/home/lcastillo98/Documents/git_projects/sim_emd_9882.mrc",
                                      "/home/lcastillo98/Documents/git_projects/emd_9882.map",
                                      segments_graph1, segments_graph2, figure1_shape, figure2_shape, './exit_fit',
                                      50)
    print("Fin")


if __name__ == '__main__':
    # main_1()
    # main_2()
    # main_3()
    main_4()
