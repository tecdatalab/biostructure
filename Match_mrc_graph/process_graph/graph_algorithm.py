import networkx as nx
import numpy as np
from operator import itemgetter

def get_total_neighbors(G, node, g_degree):
    total_neighbors = 0
    for n in G.neighbors(node):
        val_degree = 0
        for i in g_degree:
            if i[0]== n:
                val_degree = i[1]
                break
        total_neighbors += val_degree
    return total_neighbors

def similarity_function(G1, G2, g1_degree, g2_degree, node1, node2):#Que tan similares son los nodos
    node1_total_neighbors = get_total_neighbors(G1, node1, g1_degree)
    node2_total_neighbors = get_total_neighbors(G2, node2, g2_degree)
    
    return float(min(node1_total_neighbors,node2_total_neighbors))/float(max(node1_total_neighbors,node2_total_neighbors))
    
def compute_aligning_costs(G1,G2):
    g1_degree = list(G1.degree())
    g2_degree = list(G2.degree())
    
    max_degree_G1 = sorted([d for _, d in g1_degree], reverse=True)[0]
    max_degree_G2 = sorted([d for _, d in g2_degree], reverse=True)[0]
    
    C  = np.zeros((G1.number_of_nodes(), G2.number_of_nodes()))
    Z  = np.zeros((G1.number_of_nodes(), G2.number_of_nodes()))
    for i in range(G1.number_of_nodes()):
        for j in range(G2.number_of_nodes()):
            temp_center = float(g1_degree[i][1] + g2_degree[j][1])/float(max_degree_G1+max_degree_G2)
            temp_right =  float(similarity_function(G1, G2, g1_degree, g2_degree, g1_degree[i][0], g2_degree[j][0]))
            C[i][j] = (temp_center + temp_right)/2
            Z[i][j] = temp_right
    
    return C, Z

def get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, node1, node2, min_val):
    result = []
    
    g2_cylce = np.setdiff1d(list(G2.neighbors(node2)), g2_mapping)
    for n in np.setdiff1d(list(G1.neighbors(node1)), g1_mapping): #optimizar con numpy.setdiff1d
        for m in g2_cylce:
            result.append([Z[g1_dic[n]][g2_dic[m]],g1_dic[n],g2_dic[m],n,m])
    g1_check = []
    g2_check = []
    final_result = []
    total_g1_neighbors = []
    total_g2_neighbors = []
    
    result = sorted(result, key=itemgetter(0), reverse=True)
    for i in result:
        if i[0]>=min_val and i[1] not in g1_check and i[2] not in g2_check:
            final_result.append([i[1],i[2]])
            g1_check.append(i[1])
            g2_check.append(i[2])
            
        total_g1_neighbors.append(i[3])
        total_g2_neighbors.append(i[4])
    return total_g1_neighbors, total_g2_neighbors, final_result

def get_dicc_number(g_nodes):
    g_dic = {}
    for i in range(len(g_nodes)):
        g_dic[g_nodes[i]] = i
    return g_dic

def set_matrix_ij_val(matriz,i,j,val):
    matriz[i,:]=val
    matriz[:,j]=val
    return matriz

def ok_neighbors(G1, G2, g1_mapping, g2_mapping, node1, node2, equivalences_list):
    g1_check = np.intersect1d(list(G1.neighbors(node1)), g1_mapping)
    g2_check = np.intersect1d(list(G2.neighbors(node2)), g2_mapping)
    
    g2_equivalent_in_g1 = []
    for i in range(len(g2_check)): #Optimizar
        for j in range(len(equivalences_list)):
            if equivalences_list[j][1]==g2_check[i]:
                g2_equivalent_in_g1.append(equivalences_list[j][0])
                break
    
    return len(np.intersect1d(g1_check,g2_equivalent_in_g1))

def compute_neighbors_matriz(G1, G2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors, equivalences_list):
    N = np.zeros((len(total_g1_neighbors), len(total_g2_neighbors)))
    for i in range(len(total_g1_neighbors)):
        for j in range(len(total_g2_neighbors)):
            N[i][j] = ok_neighbors(G1, G2, g1_mapping, g2_mapping, total_g1_neighbors[i], total_g1_neighbors[j], equivalences_list)
    
    return N

def get_best_ij(i_list, j_list, Z, g1_dic, g2_dic, total_g1_neighbors, total_g2_neighbors):
    actual_i = i_list[0]
    actual_j = j_list[0]
    actual_val = Z[g1_dic[total_g1_neighbors[i_list[0]]]][g2_dic[total_g2_neighbors[j_list[0]]]]
    for i in range(1,len(i_list)):
        x = Z[g1_dic[total_g1_neighbors[i_list[i]]]][g2_dic[total_g2_neighbors[j_list[i]]]]
        if actual_val< x:
            actual_val = x
            actual_i = i_list[i]
            actual_j = j_list[i]
    return actual_i, actual_j

def graph_aligning(G1, G2, min_similarity_value):
    g1_nodes = list(G1.nodes)
    g2_nodes = list(G2.nodes)
    g1_dic = get_dicc_number(g1_nodes)
    g2_dic = get_dicc_number(g2_nodes)
        
    C, Z = compute_aligning_costs(G1,G2)
    result = []
    
    g1_mapping = []
    g2_mapping = []
    total_g1_neighbors = []
    total_g2_neighbors = []
    max_similarity_value = np.amax(C)
    
    while max_similarity_value>= min_similarity_value:
        i,j = np.where(C == max_similarity_value)
        '''Cambiar esto por escogencia mas cerca'''
        result.append([g1_nodes[i[0]],g2_nodes[j[0]]])
        g1_mapping.append(g1_nodes[i[0]])
        g2_mapping.append(g2_nodes[j[0]])
        C = set_matrix_ij_val(C,i,j,-1)
        Z = set_matrix_ij_val(Z,i,j,-1)
         
        g1_neighbors ,g2_neighbors, queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, g1_nodes[i[0]], g2_nodes[j[0]],min_similarity_value)
        total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
        total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
        while queue_neighbors:
            temp = queue_neighbors.pop(0)
            i,j = temp[0], temp[1]
            if Z[i][j]==-1:
                continue
            result.append([g1_nodes[i],g2_nodes[j]])
            g1_mapping.append(g1_nodes[i])
            g2_mapping.append(g2_nodes[j])
            C = set_matrix_ij_val(C,i,j,-1)
            Z = set_matrix_ij_val(Z,i,j,-1)
            
            g1_neighbors ,g2_neighbors, new_queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, g1_nodes[i], g2_nodes[j],min_similarity_value)
            total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
            total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
            queue_neighbors+= new_queue_neighbors
        
        max_similarity_value = np.amax(C)
    '''Arreglar el problema de los valores externos'''
    
    total_g1_neighbors = np.setdiff1d(total_g1_neighbors, g1_mapping)
    total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping) 
    if len(total_g1_neighbors)> 0 and len(total_g2_neighbors)> 0: # Todo esta ok si no entra
        N = compute_neighbors_matriz(G1, G2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors, result)
        max_similarity_value = np.amax(N)
        
        while max_similarity_value>0:
            i,j = np.where(N == max_similarity_value)
            i,j = get_best_ij(i, j, Z, g1_dic, g2_dic, total_g1_neighbors, total_g2_neighbors)
            
            result.append([total_g1_neighbors[i],total_g2_neighbors[j]])
            g1_mapping.append(total_g1_neighbors[i])
            g2_mapping.append(total_g2_neighbors[j])
            C = set_matrix_ij_val(C,g1_dic[total_g1_neighbors[i]],g2_dic[total_g2_neighbors[j]],-1)
            Z = set_matrix_ij_val(Z,g1_dic[total_g1_neighbors[i]],g2_dic[total_g2_neighbors[j]],-1)
            
            g1_neighbors ,g2_neighbors, new_queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, total_g1_neighbors[i], total_g2_neighbors[j],0)
            total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
            total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
            
            total_g1_neighbors = np.setdiff1d(total_g1_neighbors, g1_mapping)
            total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping)
    return result



