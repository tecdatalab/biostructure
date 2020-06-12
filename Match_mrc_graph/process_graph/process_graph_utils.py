'''
Created on Jun 9, 2020

@author: luis98
'''
from sklearn.neighbors import KDTree
from process_graph.process_segment_faces import get_n_points_cube
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

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

