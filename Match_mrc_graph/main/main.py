from process_mrc.miscellaneou import get_mrc_segments, \
    get_mrc_synthetic_segments_pdb
from process_graph.process_graph_utils import generate_graph,\
    draw_graph_similarity, draw_graph_similarity_same_image,\
    get_similarity_complete_cloud
from process_graph.graph_algorithm import graph_aligning
from fit.fit_map import fit_map_in_map

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
    result = fit_map_in_map('/home/luis98/Documents/pruebas/175d_A.mrc', '/home/luis98/Documents/pruebas/175d.mrc')
    result.print_data()

if __name__ == '__main__':
    #main()
    main_1()

