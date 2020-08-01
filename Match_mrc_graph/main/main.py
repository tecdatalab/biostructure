import shutil

from process_mrc.generate import get_mrc_segments, \
    get_mrc_synthetic_segments_pdb
from process_graph.process_graph_utils import generate_graph, \
    draw_graph_similarity, draw_graph_similarity_same_image, \
    get_similarity_complete_cloud
from process_graph.graph_algorithm import graph_aligning
from fit.fit_map import fit_map_in_map
from pdb_mrc.pdb_transform import download_pdb, get_chains, pdb_to_mrc_chains
import os

from process_mrc.miscellaneous import get_center_point
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


def main_1():
    path = 'temp_data'
    try:
        os.mkdir(path)
        download_pdb('175d', './temp_data/175d.pdb')
        download_pdb('6m03', './temp_data/6m03.pdb')
    except:
        pass

    # Creacion de los mapas
    chains = get_chains('./temp_data/175d.pdb')
    pdb_to_mrc_chains(True, False, 5.0, './temp_data/175d.pdb', './temp_data', chains, 1)

    chains = get_chains('./temp_data/6m03.pdb')
    pdb_to_mrc_chains(True, False, 5.0, './temp_data/6m03.pdb', './temp_data', chains, 1)

    # Generacion de grafos a partir de la segmentacion de los mapas
    # segments_graph1 = get_mrc_synthetic_segments_pdb("./temp_data/175d", 7)
    # segments_graph2 = get_mrc_synthetic_segments_pdb("./temp_data/6m03", 7)

    segments_graph1 = get_mrc_segments("./temp_data/175d/175d.mrc", 7, 3, 1)
    segments_graph2 = get_mrc_segments("./temp_data/6m03/6m03.mrc", 7, 3, 1)

    print("Segments 1 shape", segments_graph1[0].mask.shape)
    print("Segments 2 shape", segments_graph2[0].mask.shape)

    graph1 = generate_graph(segments_graph1, 50, 0, 6, 1)
    graph2 = generate_graph(segments_graph2, 50, 0, 6, 1)

    # Match de grafos
    result = graph_aligning(graph1, graph2, 1, False)
    print("Result match graph", result)

    graph1_match_index = get_element_list(0, result)
    graph2_match_index = get_element_list(1, result)

    # Puntos centricos
    # print(graph1_match_index)
    # print(graph2_match_index)

    center_point1 = get_center_point(graph1_match_index, segments_graph1, 0)
    center_point2 = get_center_point(graph2_match_index, segments_graph2, 0)

    print("Center point 1:", center_point1)
    print("Center point 2:", center_point2)

    # Fit map in map
    result = fit_map_in_map('./temp_data/6m03/6m03.mrc', './temp_data/175d/175d.mrc', './exit_fit')
    result.print_data()


if __name__ == '__main__':
    # main()
    main_1()
