import networkx as nx
import numpy as np
from operator import itemgetter
from matplotlib.pyplot import axis
from sklearn.metrics import mean_squared_error

def get_total_neighbors_of_neighbors(G, node, g_degree): #Obtiene la cantidad de vecinos para todos los vecinos de un nodo dado
    total_neighbors = 0
    for n in G.neighbors(node):
        val_degree = 0
        for i in g_degree:
            if i[0]== n:
                val_degree = i[1]
                break
        total_neighbors += val_degree
    return total_neighbors

def setdiff1d_multidimentional(list_A, list_B):
    a1_rows = set(map(tuple, list_A))
    a2_rows = set(map(tuple, list_B))
    result = [list(elem) for elem in list(a1_rows.difference(a2_rows))]
    return result

def structural_similarity_node(G1, G2, g1_degree, g2_degree, node1, node2):#Que tan similares son los nodos
    node1_total_neighbors = get_total_neighbors_of_neighbors(G1, node1, g1_degree)
    node2_total_neighbors = get_total_neighbors_of_neighbors(G2, node2, g2_degree)
    
    similarity = 1.0 - float(min(node1_total_neighbors,node2_total_neighbors))/float(max(node1_total_neighbors,node2_total_neighbors))
    return similarity #Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen

def similarity_node(g1_nodes_with_data, g2_nodes_with_data, node1, node2):#Que tan similares son los nodos mediantre el error minimo cuadrado
    mse = mean_squared_error(g1_nodes_with_data[node1]['zd_descriptors'], g2_nodes_with_data[node2]['zd_descriptors'])
    return mse #Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen
    
def compute_aligning_costs(G1, G2, g1_nodes_with_data, g2_nodes_with_data):
    g1_degree = list(G1.degree()) #Cantidad de aristas por nodo
    g2_degree = list(G2.degree())
    
    max_degree_G1 = sorted([d for _, d in g1_degree], reverse=True)[0] #Creo que se puede mejorar esto
    max_degree_G2 = sorted([d for _, d in g2_degree], reverse=True)[0] #Creo que se puede mejorar esto
    
    C  = np.zeros((G1.number_of_nodes(), G2.number_of_nodes())) #Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
    Z  = np.zeros((G1.number_of_nodes(), G2.number_of_nodes()))
    
    for i in range(G1.number_of_nodes()):
        for j in range(G2.number_of_nodes()): # Para todos los nodos calcular su similidad (Z) y potencial de convertirse en raiz (C)
            
            nodes_str_simi = 1.0 - float(g1_degree[i][1] + g2_degree[j][1])/float(max_degree_G1+max_degree_G2)
            total_str_simi =  float(structural_similarity_node(G1, G2, g1_degree, g2_degree, g1_degree[i][0], g2_degree[j][0]))
            sim_node = similarity_node(g1_nodes_with_data, g2_nodes_with_data, g1_degree[i][0], g2_degree[j][0])
            #print(similarity_node(G1, G2, g1_degree[i][0], g2_degree[j][0]))
            C[i][j] = nodes_str_simi + total_str_simi + sim_node
            #Z[i][j] = total_str_simi
            Z[i][j] = sim_node
    
    return C, Z

def get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, node1, node2, min_val):
    result = []
    
    g2_cylce = np.setdiff1d(list(G2.neighbors(node2)), g2_mapping) #Se optienen los vecinos que no han sido mapeados ya
    for n in np.setdiff1d(list(G1.neighbors(node1)), g1_mapping): 
        for m in g2_cylce:
            result.append([Z[g1_dic[n]][g2_dic[m]],g1_dic[n],g2_dic[m],n,m]) #A lista de resultado se le agregan los valores de: [Similidad, index n, index m, n value, m value]
            
    g1_check = [] #Lista para revisar que un elemento no ha sido ya mapedo con otro nodo
    g2_check = []
    final_result = [] #Lista final de parejas
    total_g1_neighbors = []#Vecinos de los nodos vistos
    total_g2_neighbors = []
    
    result = sorted(result, key=itemgetter(0), reverse=False)# Se orden de menor a mayor los elementos a agregar
    for i in result:
        if i[0]<=min_val and i[1] not in g1_check and i[2] not in g2_check:
            final_result.append([i[1],i[2]])
            g1_check.append(i[1])
            g2_check.append(i[2])
            
        total_g1_neighbors.append(i[3])
        total_g2_neighbors.append(i[4])
    return total_g1_neighbors, total_g2_neighbors, final_result

def get_dicc_number(g_nodes): #Dado un alista de nodos, se crea un diccionario llave -> key [0,1...]
    g_dic = {}
    for i in range(len(g_nodes)):
        g_dic[g_nodes[i]] = i
    return g_dic

def set_matrix_ij_val(matriz,i,j,val): #Se setea una fila y columna dada con un valor dado
    matriz[i,:]=val
    matriz[:,j]=val
    return matriz

def ok_neighbors(G1, G2, g1_mapping, g2_mapping, node1, node2, equivalences_list): #Dado dos nodos se retorna la cantidad de nodos compartidos, mapeados
    g1_check = np.intersect1d(list(G1.neighbors(node1)), g1_mapping)
    g2_check = np.intersect1d(list(G2.neighbors(node2)), g2_mapping)
    
    g2_equivalent_in_g1 = []
    for i in range(len(g2_check)): #Optimizar
        for j in range(len(equivalences_list)):
            if equivalences_list[j][1]==g2_check[i]:
                g2_equivalent_in_g1.append(equivalences_list[j][0])
                break
    
    return len(np.intersect1d(g1_check,g2_equivalent_in_g1))

def compute_neighbors_matriz(G1, G2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors, equivalences_list): #Calculo de matriz de nodos compartidos
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

def graph_aligning(G1, G2, min_similarity_value, repair_external_values = True):
    g1_nodes = list(G1.nodes) #Obtener lista de nodos
    g2_nodes = list(G2.nodes) 
    g1_nodes_with_data = dict(G1.nodes(data=True)) #Obtener lista de nodos
    g2_nodes_with_data = dict(G2.nodes(data=True)) 
    g1_dic = get_dicc_number(g1_nodes) #A cada nodo se le asigna un indice, con el cual estan almacendos en la matriz
    g2_dic = get_dicc_number(g2_nodes)
        
    C, Z = compute_aligning_costs(G1, G2, g1_nodes_with_data, g2_nodes_with_data)#Computar matriz con la capacidad de convertirse en raiz(C) o similaridad (Z)
    result = []
    
    g1_mapping = [] #Se almacenan todos los nodos que fueron mapeados (con su valor real)
    g2_mapping = []
    total_g1_neighbors = [] #Se alamacenan todos los vecinos externos 
    total_g2_neighbors = []
    most_similarity_value = np.amin(C)
    while most_similarity_value<= min_similarity_value:
        i,j = np.where(C == most_similarity_value)
        '''Cambiar esto por escogencia probilistica'''
        result.append([g1_nodes[i[0]],g2_nodes[j[0]]]) #Se agrega el nuevo match
        g1_mapping.append(g1_nodes[i[0]]) #Se agrega el nuevo match los mapedos
        g2_mapping.append(g2_nodes[j[0]])
        C = set_matrix_ij_val(C,i,j,np.Inf) #Se llenan con infinito las filas y columnas que no pueden ser utilizadas
        Z = set_matrix_ij_val(Z,i,j,np.Inf)
         
        g1_neighbors ,g2_neighbors, queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, g1_nodes[i[0]], g2_nodes[j[0]],min_similarity_value) # Obtienen los vecinos de los nodos mapeados que son mayores a un umbral
        total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist() #Se agrega a la lista total de vecinos, los vecinos que no han sido agregados
        total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
        while queue_neighbors: #Crecemos en ancho para cada un de los vecinos
            temp = queue_neighbors.pop(0)
            i,j = temp[0], temp[1]
            if Z[i][j]==np.Inf:
                continue
            result.append([g1_nodes[i],g2_nodes[j]])
            g1_mapping.append(g1_nodes[i])
            g2_mapping.append(g2_nodes[j])
            C = set_matrix_ij_val(C,i,j,np.Inf)
            Z = set_matrix_ij_val(Z,i,j,np.Inf)
            
            g1_neighbors ,g2_neighbors, new_queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, g1_nodes[i], g2_nodes[j],min_similarity_value)
            total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
            total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
            queue_neighbors += setdiff1d_multidimentional(new_queue_neighbors, queue_neighbors) 
        
        most_similarity_value = np.amin(C)
    
    if not repair_external_values:return result
    '''Arreglar el problema de los valores externos'''
    total_g1_neighbors = np.setdiff1d(total_g1_neighbors, g1_mapping) # Se limpian los vecinos con los elementos que ya fueron mapeados
    total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping) 
    if len(total_g1_neighbors)> 0 and len(total_g2_neighbors)> 0: # Si existen vecinos libres para ambos lados
        N = compute_neighbors_matriz(G1, G2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors, result)
        max_similarity_value = np.amax(N) #Se obtiene el elemento con mas nodos compartidos
        
        while max_similarity_value>0:
            i,j = np.where(N == max_similarity_value)
            i,j = get_best_ij(i, j, Z, g1_dic, g2_dic, total_g1_neighbors, total_g2_neighbors) #De todos los posibles valores a hacer match, hacer match con el mejor
            
            result.append([total_g1_neighbors[i],total_g2_neighbors[j]])
            g1_mapping.append(total_g1_neighbors[i])
            g2_mapping.append(total_g2_neighbors[j])
            C = set_matrix_ij_val(C,g1_dic[total_g1_neighbors[i]],g2_dic[total_g2_neighbors[j]],np.Inf)
            Z = set_matrix_ij_val(Z,g1_dic[total_g1_neighbors[i]],g2_dic[total_g2_neighbors[j]],np.Inf)
            
            g1_neighbors ,g2_neighbors, new_queue_neighbors = get_promising_neighbors(G1, G2, g1_dic, g2_dic, g1_mapping, g2_mapping, Z, total_g1_neighbors[i], total_g2_neighbors[j],np.Inf)
            total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
            total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
            
            total_g1_neighbors = np.setdiff1d(total_g1_neighbors, g1_mapping)
            total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping)
    return result



