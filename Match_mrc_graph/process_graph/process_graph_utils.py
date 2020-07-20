'''
Created on Jun 9, 2020

@author: luis98
'''
from sklearn.neighbors import KDTree
from process_graph.process_segment_faces import get_n_points_cube
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def generate_graph(segments, n_points_face, filter_value, max_distance, min_points):
    face_points = {}
    kd_trees = {}
    g_result=nx.Graph()
    
    for i in segments:
        face_points[i.id_segment] = get_n_points_cube(i.mask, n_points_face, filter_value)
        
    for i in segments:
        kd_trees[i.id_segment] = KDTree(face_points[i.id_segment])
        
    check = []
    for i in segments:
        for j in segments:
            if (i.id_segment != j.id_segment) and ([i.id_segment,j.id_segment] not in check):
                #check code
                ok_values = 0
                for point in face_points[i.id_segment]:
                    dist, _ind = kd_trees[j.id_segment].query([point], k=1) 
                    if dist[0]<=max_distance:
                        ok_values+=1
                        if ok_values>=min_points:
                            ### Creation of graph
                            g_result.add_node(i.id_segment, zd_descriptors=i.zd_descriptors)
                            g_result.add_node(j.id_segment, zd_descriptors=j.zd_descriptors)
                            g_result.add_edge(i.id_segment, j.id_segment)
                            check.append([j.id_segment,i.id_segment])
                            break
                check.append([i.id_segment,j.id_segment])
            
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
    nx.draw_networkx(graph1, ax=ax[2], with_labels=True, node_color=graph1_match)
    ax[2].set_axis_off()
    
    ax[3].set_title('Graph 1 match')
    nx.draw_networkx(graph2, ax=ax[3], with_labels=True, node_color=graph2_match)
    ax[3].set_axis_off()

    plt.show()

def get_similarity_complete_cloud(graph1, graph2, similary):
    list_points_graph1 = []
    list_points_graph2 = []
    for i in similary:
        for k in graph1:
            if k.id_segment == i[0]:
                x,y,z = np.where(k.mask>0)
                for w in range(len(x)):
                    list_points_graph1.append([x[w],y[w],z[w]])
                print(len(x))
                break
        for k in graph2:
            if k.id_segment == i[1]:
                x,y,z = np.where(k.mask>0)
                for w in range(len(x)):
                    list_points_graph2.append([x[w],y[w],z[w]])
                break
    
    
    
    return (np.array(list_points_graph1),np.array(list_points_graph2))





