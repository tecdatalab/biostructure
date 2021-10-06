'''
Created on Jun 9, 2020

@author: luis98
'''
import itertools

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from sklearn.neighbors import KDTree
from process_graph.process_segment_faces import get_n_points_cube
import threading

from process_mrc.miscellaneous import get_sum_xyz_can


def gen_thread_points(id, mask, n_points_face, filter_value, semaphore, dic_result):
  result = get_n_points_cube(mask, n_points_face, filter_value)
  semaphore.acquire()
  dic_result[id] = result
  semaphore.release()


def gen_thread_KDTree(id, data, semaphore, dic_result):
  result = KDTree(data)
  semaphore.acquire()
  dic_result[id] = result
  semaphore.release()


def generate_graph(segments, n_points_face, filter_value, max_distance, min_points):
  segments_id = []
  for i in segments:
    if i.id_segment in segments_id:
      raise TypeError("Error in segments id, some are duplicate")
    else:
      segments_id.append(i.id_segment)

  face_points = {}
  kd_trees = {}
  g_result = nx.Graph()

  sem = threading.Semaphore()
  thread_list = []
  for i in segments:
    t = threading.Thread(target=gen_thread_points, args=[i.id_segment, i.mask, n_points_face, filter_value,
                                                         sem, face_points])

    thread_list.append(t)
    t.start()

  for t in thread_list:
    t.join()

  thread_list = []
  for i in segments:
    t = threading.Thread(target=gen_thread_KDTree, args=[i.id_segment, face_points[i.id_segment], sem,
                                                         kd_trees])
    thread_list.append(t)
    t.start()

  for t in thread_list:
    t.join()

    # Add nodes in gaphs
  for i in segments:
    g_result.add_node(i.id_segment,
                      zd_descriptors=i.zd_descriptors,
                      cube_xyz_can=get_sum_xyz_can(i),
                      volume=i.volume,
                      textSimulatePDB=i.textSimulatePDB)

  for item in itertools.combinations(np.arange(len(segments)), 2):
    if item[0] == item[1]:
      continue

    i = segments[item[0]]
    j = segments[item[1]]

    # check code
    dist, _ind = kd_trees[j.id_segment].query(face_points[i.id_segment], k=1)
    chek_val = np.max(np.sort(np.squeeze(dist))[:min(len(dist), min_points)])

    if chek_val <= max_distance:
      g_result.add_edge(i.id_segment, j.id_segment)

  return g_result


def draw_graph_similarity(graph1, graph2, result):
  matplotlib.use('TKAgg')
  graph1_nodes = []
  graph2_nodes = []

  for i in result:
    graph1_nodes.append(i[0])
    graph2_nodes.append(i[1])

  graph1_match = []
  graph2_match = []

  for i in list(graph1.nodes()):
    if i in graph1_nodes:
      graph1_match.append('y')
    else:
      graph1_match.append('r')

  for i in list(graph2.nodes()):
    if i in graph2_nodes:
      graph2_match.append('y')
    else:
      graph2_match.append('b')

  plt.figure(1)
  plt.suptitle('Graph 0', fontsize=16)
  nx.draw(graph1, with_labels=True, node_color='r')
  plt.figure(2)
  plt.suptitle('Graph 1', fontsize=16)
  nx.draw(graph2, with_labels=True, node_color='b')
  plt.figure(3)
  plt.suptitle('Graph 0 match', fontsize=16)
  nx.draw(graph1, with_labels=True, node_color=graph1_match)
  plt.figure(4)
  plt.suptitle('Graph 1 match', fontsize=16)
  nx.draw(graph1, with_labels=True, node_color=graph2_match)
  plt.show()


def draw_graph_similarity_same_image(graph1, graph2, result):
  matplotlib.use('TKAgg')

  graph1_nodes = []
  graph2_nodes = []

  graph1_match = []
  graph2_match = []

  for i in result:
    graph1_nodes.append(i[0])
    graph2_nodes.append(i[1])

  for i in list(graph1.nodes()):
    if i in graph1_nodes:
      graph1_match.append('y')
    else:
      graph1_match.append('r')

  for i in list(graph2.nodes()):
    if i in graph2_nodes:
      graph2_match.append('y')
    else:
      graph2_match.append('b')

  fig, axes = plt.subplots(nrows=2, ncols=2)
  ax = axes.flatten()

  ax[0].set_title('Graph 0')
  nx.draw_networkx(graph1, ax=ax[0], with_labels=True, node_color='r')
  ax[0].set_axis_off()

  ax[1].set_title('Graph 1')
  nx.draw_networkx(graph2, ax=ax[1], with_labels=True, node_color='b')
  ax[1].set_axis_off()

  ax[2].set_title('Graph 0 match')
  mapping = {}
  for i in result:
    mapping[i[0]] = str(i[0]) + ":" + str(i[1])
  graph1_cpy = nx.relabel_nodes(graph1, mapping)
  nx.draw_networkx(graph1_cpy, ax=ax[2], with_labels=True, node_color=graph1_match)
  ax[2].set_axis_off()

  ax[3].set_title('Graph 1 match')
  mapping = {}
  for i in result:
    mapping[i[1]] = str(i[1]) + ":" + str(i[0])
  graph2_cpy = nx.relabel_nodes(graph2, mapping)
  nx.draw_networkx(graph2_cpy, ax=ax[3], with_labels=True, node_color=graph2_match)
  ax[3].set_axis_off()

  plt.show()


def get_similarity_complete_cloud(graph1, graph2, similary):
  list_points_graph1 = []
  list_points_graph2 = []
  for i in similary:
    for k in graph1:
      if k.id_segment == i[0]:
        x, y, z = np.where(k.mask > 0)
        for w in range(len(x)):
          list_points_graph1.append([x[w], y[w], z[w]])
        print(len(x))
        break
    for k in graph2:
      if k.id_segment == i[1]:
        x, y, z = np.where(k.mask > 0)
        for w in range(len(x)):
          list_points_graph2.append([x[w], y[w], z[w]])
        break

  return (np.array(list_points_graph1), np.array(list_points_graph2))
