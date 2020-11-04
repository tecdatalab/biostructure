from operator import itemgetter

import numpy as np
from sklearn.metrics import mean_squared_error


def get_total_neighbors_of_neighbors(graph, node,
                                     g_degree):  # Obtiene la cantidad de vecinos para todos los vecinos de un nodo dado
  total_neighbors = 0
  for n in graph.neighbors(node):
    val_degree = 0
    for i in g_degree:
      if i[0] == n:
        val_degree = i[1]
        break
    total_neighbors += val_degree
  return total_neighbors


def set_diff1d_multidimensional(list_a, list_b):
  a1_rows = set(map(tuple, list_a))
  a2_rows = set(map(tuple, list_b))
  result = [list(elem) for elem in list(a1_rows.difference(a2_rows))]
  return result


def structural_similarity_node(graph1, graph2, g1_degree, g2_degree, node1, node2):  # Que tan similares son los nodos
  node1_total_neighbors = get_total_neighbors_of_neighbors(graph1, node1, g1_degree)
  node2_total_neighbors = get_total_neighbors_of_neighbors(graph2, node2, g2_degree)

  if float(min(node1_total_neighbors, node2_total_neighbors)) == float(max(node1_total_neighbors, node2_total_neighbors)):
    sim_error = 1.0
  else:
    if max(node1_total_neighbors, node2_total_neighbors) == 0:
      sim_error = 0
    else:
      sim_error = \
      float(min(node1_total_neighbors, node2_total_neighbors)) / float(max(node1_total_neighbors, node2_total_neighbors))

  similarity = 1.0 - sim_error
  return similarity  # Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen


def similarity_node(g1_nodes_with_data, g2_nodes_with_data, node1,
                    node2):  # Que tan similares son los nodos mediantre el error minimo cuadrado
  mse = mean_squared_error(g1_nodes_with_data[node1]['zd_descriptors'], g2_nodes_with_data[node2]['zd_descriptors'])
  return mse  # Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen


def compute_aligning_costs(graph1, graph2, g1_nodes_with_data, g2_nodes_with_data):
  g1_degree = list(graph1.degree())  # Cantidad de aristas por nodo
  g2_degree = list(graph2.degree())

  # max_degree_g1 = sorted([d for _, d in g1_degree], reverse=True)[0]  # Creo que se puede mejorar esto
  max_degree_g1 = np.max(np.array(g1_degree).reshape(2, len(g1_degree))[1])
  # max_degree_g2 = sorted([d for _, d in g2_degree], reverse=True)[0]  # Creo que se puede mejorar esto
  max_degree_g2 = np.max(np.array(g2_degree).reshape(2, len(g2_degree))[1])

  C = np.zeros((graph1.number_of_nodes(),
                graph2.number_of_nodes()))  # Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
  Z = np.zeros((graph1.number_of_nodes(), graph2.number_of_nodes()))

  for i in range(graph1.number_of_nodes()):
    for j in range(
      graph2.number_of_nodes()):  # Para todos los nodos calcular su similidad (Z) y potencial de convertirse en raiz (C)

      if float(g1_degree[i][1] + g2_degree[j][1]) == float(max_degree_g1 + max_degree_g2):
        sim_val = 1.0
      else:
        if max_degree_g1 + max_degree_g2 == 0:
          sim_val = 0
        else:
          sim_val = float(g1_degree[i][1] + g2_degree[j][1]) / float(max_degree_g1 + max_degree_g2)
      nodes_slr_sim = 1.0 - sim_val
      total_str_sim = \
        float(structural_similarity_node(graph1, graph2, g1_degree, g2_degree, g1_degree[i][0], g2_degree[j][0]))
      sim_node = similarity_node(g1_nodes_with_data, g2_nodes_with_data, g1_degree[i][0], g2_degree[j][0])
      # print(similarity_node(G1, G2, g1_degree[i][0], g2_degree[j][0]))
      C[i][j] = nodes_slr_sim + total_str_sim + sim_node
      # Z[i][j] = total_str_sim
      Z[i][j] = sim_node

  return C, Z


def get_promising_neighbors(graph1, graph2, g1_dic, g2_dic, g1_mapping, g2_mapping, matrix_z, node1, node2, min_val):
  result = []

  g2_cycle = np.setdiff1d(list(graph2.neighbors(node2)),
                          g2_mapping)  # Se optienen los vecinos que no han sido mapeados ya
  for n in np.setdiff1d(list(graph1.neighbors(node1)), g1_mapping):
    for m in g2_cycle:
      result.append([matrix_z[g1_dic[n]][g2_dic[m]], g1_dic[n], g2_dic[m], n,
                     m])  # A lista de resultado se le agregan los valores de: [Similidad, index n, index m, n value, m value]

  g1_check = []  # Lista para revisar que un elemento no ha sido ya mapedo con otro nodo
  g2_check = []
  final_result = []  # Lista final de parejas
  total_g1_neighbors = []  # Vecinos de los nodos vistos
  total_g2_neighbors = []

  result = sorted(result, key=itemgetter(0), reverse=False)  # Se orden de menor a mayor los elementos a agregar
  for i in result:
    if i[0] <= min_val and i[1] not in g1_check and i[2] not in g2_check:
      final_result.append([i[1], i[2]])
      g1_check.append(i[1])
      g2_check.append(i[2])

    total_g1_neighbors.append(i[3])
    total_g2_neighbors.append(i[4])
  return total_g1_neighbors, total_g2_neighbors, final_result


def get_dic_number(g_nodes):  # Dado un alista de nodos, se crea un diccionario llave -> key [0,1...]
  g_dic = {}
  for i in range(len(g_nodes)):
    g_dic[g_nodes[i]] = i
  return g_dic


def set_matrix_ij_val(matrix, i, j, val):  # Se setea una fila y columna dada con un valor dado
  matrix[i, :] = val
  matrix[:, j] = val
  return matrix


def ok_neighbors(graph1, graph2, g1_mapping, g2_mapping, node1, node2,
                 equivalences_list):  # Dado dos nodos se retorna la cantidad de nodos compartidos, mapeados
  g1_check = np.intersect1d(list(graph1.neighbors(node1)), g1_mapping)
  g2_check = np.intersect1d(list(graph2.neighbors(node2)), g2_mapping)

  g2_equivalent_in_g1 = []
  for i in range(len(g2_check)):  # Optimizar
    for j in range(len(equivalences_list)):
      if equivalences_list[j][1] == g2_check[i]:
        g2_equivalent_in_g1.append(equivalences_list[j][0])
        break

  return len(np.intersect1d(g1_check, g2_equivalent_in_g1))


def compute_neighbors_matrix(graph1, graph2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors,
                             equivalences_list):  # Calculo de matriz de nodos compartidos
  N = np.zeros((len(total_g1_neighbors), len(total_g2_neighbors)))
  for i in range(len(total_g1_neighbors)):
    for j in range(len(total_g2_neighbors)):
      N[i][j] = ok_neighbors(graph1, graph2, g1_mapping, g2_mapping, total_g1_neighbors[i], total_g1_neighbors[j],
                             equivalences_list)

  return N


def get_best_ij(i_list, j_list, matrix_z, g1_dic, g2_dic, total_g1_neighbors, total_g2_neighbors):
  actual_i = i_list[0]
  actual_j = j_list[0]
  actual_val = matrix_z[g1_dic[total_g1_neighbors[i_list[0]]]][g2_dic[total_g2_neighbors[j_list[0]]]]
  for i in range(1, len(i_list)):
    x = matrix_z[g1_dic[total_g1_neighbors[i_list[i]]]][g2_dic[total_g2_neighbors[j_list[i]]]]
    if actual_val < x:
      actual_val = x
      actual_i = i_list[i]
      actual_j = j_list[i]
  return actual_i, actual_j


def graph_aligning(graph1, graph2, min_similarity_value, repair_external_values=True):
  g1_nodes = list(graph1.nodes)  # Obtener lista de nodos
  g2_nodes = list(graph2.nodes)
  g1_nodes_with_data = dict(graph1.nodes(data=True))  # Obtener lista de nodos
  g2_nodes_with_data = dict(graph2.nodes(data=True))
  g1_dic = get_dic_number(g1_nodes)  # A cada nodo se le asigna un indice, con el cual estan almacendos en la matriz
  g2_dic = get_dic_number(g2_nodes)

  C, Z = compute_aligning_costs(graph1, graph2, g1_nodes_with_data,
                                g2_nodes_with_data)  # Computar matriz con la capacidad de convertirse en raiz(C) o similaridad (Z)
  result = []

  g1_mapping = []  # Se almacenan todos los nodos que fueron mapeados (con su valor real)
  g2_mapping = []
  total_g1_neighbors = []  # Se alamacenan todos los vecinos externos
  total_g2_neighbors = []
  most_similarity_value = np.amin(C)
  while most_similarity_value <= min_similarity_value:
    i, j = np.where(C == most_similarity_value)
    '''Cambiar esto por escogencia probilistica'''
    result.append([g1_nodes[i[0]], g2_nodes[j[0]]])  # Se agrega el nuevo match
    g1_mapping.append(g1_nodes[i[0]])  # Se agrega el nuevo match los mapedos
    g2_mapping.append(g2_nodes[j[0]])
    C = set_matrix_ij_val(C, i, j,
                          np.Inf)  # Se llenan con infinito las filas y columnas que no pueden ser utilizadas
    Z = set_matrix_ij_val(Z, i, j, np.Inf)

    g1_neighbors, g2_neighbors, queue_neighbors = get_promising_neighbors(graph1, graph2, g1_dic, g2_dic, g1_mapping,
                                                                          g2_mapping, Z, g1_nodes[i[0]],
                                                                          g2_nodes[j[0]],
                                                                          min_similarity_value)  # Obtienen los vecinos de los nodos mapeados que son mayores a un umbral
    total_g1_neighbors += np.setdiff1d(g1_neighbors,
                                       total_g1_neighbors).tolist()  # Se agrega a la lista total de vecinos, los vecinos que no han sido agregados
    total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
    while queue_neighbors:  # Crecemos en ancho para cada un de los vecinos
      temp = queue_neighbors.pop(0)
      i, j = temp[0], temp[1]
      if Z[i][j] == np.Inf:
        continue
      result.append([g1_nodes[i], g2_nodes[j]])
      g1_mapping.append(g1_nodes[i])
      g2_mapping.append(g2_nodes[j])
      C = set_matrix_ij_val(C, i, j, np.Inf)
      Z = set_matrix_ij_val(Z, i, j, np.Inf)

      g1_neighbors, g2_neighbors, new_queue_neighbors = get_promising_neighbors(graph1, graph2, g1_dic, g2_dic,
                                                                                g1_mapping, g2_mapping, Z,
                                                                                g1_nodes[i], g2_nodes[j],
                                                                                min_similarity_value)
      total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
      total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()
      queue_neighbors += set_diff1d_multidimensional(new_queue_neighbors, queue_neighbors)

    most_similarity_value = np.amin(C)

  if not repair_external_values: return result
  '''Arreglar el problema de los valores externos'''
  total_g1_neighbors = np.setdiff1d(total_g1_neighbors,
                                    g1_mapping)  # Se limpian los vecinos con los elementos que ya fueron mapeados
  total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping)
  if len(total_g1_neighbors) > 0 and len(total_g2_neighbors) > 0:  # Si existen vecinos libres para ambos lados
    N = compute_neighbors_matrix(graph1, graph2, g1_mapping, g2_mapping, total_g1_neighbors, total_g2_neighbors, result)
    max_similarity_value = np.amax(N)  # Se obtiene el elemento con mas nodos compartidos

    while max_similarity_value > 0:
      i, j = np.where(N == max_similarity_value)
      i, j = get_best_ij(i, j, Z, g1_dic, g2_dic, total_g1_neighbors,
                         total_g2_neighbors)  # De todos los posibles valores a hacer match, hacer match con el mejor

      result.append([total_g1_neighbors[i], total_g2_neighbors[j]])
      g1_mapping.append(total_g1_neighbors[i])
      g2_mapping.append(total_g2_neighbors[j])
      C = set_matrix_ij_val(C, g1_dic[total_g1_neighbors[i]], g2_dic[total_g2_neighbors[j]], np.Inf)
      Z = set_matrix_ij_val(Z, g1_dic[total_g1_neighbors[i]], g2_dic[total_g2_neighbors[j]], np.Inf)

      g1_neighbors, g2_neighbors, new_queue_neighbors = get_promising_neighbors(graph1, graph2, g1_dic, g2_dic,
                                                                                g1_mapping, g2_mapping, Z,
                                                                                total_g1_neighbors[i],
                                                                                total_g2_neighbors[j], np.Inf)
      total_g1_neighbors += np.setdiff1d(g1_neighbors, total_g1_neighbors).tolist()
      total_g2_neighbors += np.setdiff1d(g2_neighbors, total_g2_neighbors).tolist()

      total_g1_neighbors = np.setdiff1d(total_g1_neighbors, g1_mapping)
      total_g2_neighbors = np.setdiff1d(total_g2_neighbors, g2_mapping)
  return result
