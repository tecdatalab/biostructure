from operator import itemgetter

import numpy as np
import random
from sklearn.metrics import mean_squared_error

from general_utils.graph_utils import get_node_by_position, get_node_name_by_pos, get_node_index
from general_utils.list_utils import setdiff2d, intersect2d


# Obtiene la cantidad de vecinos para todos los vecinos de un nodo dado
def cal_degree_non(graph, node):
  total_neighbors = 0
  for n in graph.neighbors(node):
    total_neighbors += graph.degree[n]
  return total_neighbors


# Que tan similares son los nodos
def structural_similarity_node(graph1, graph2, node1, node2):
  node1_total_neighbors = cal_degree_non(graph1, node1)
  node2_total_neighbors = cal_degree_non(graph2, node2)

  if node1_total_neighbors == node2_total_neighbors:
    sim_error = 0.0
  else:
    sim_error = 1.0 - \
                float(min(node1_total_neighbors, node2_total_neighbors)) / \
                float(max(node1_total_neighbors, node2_total_neighbors))

  return sim_error


# Que tan similares son los nodos mediantre el error minimo cuadrado
def similarity_node(node1, node2):
  mse = mean_squared_error(node1['zd_descriptors'], node2['zd_descriptors'])
  return mse  # Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen


def calculate_C(graph1, graph2):
  max_degree_g1 = np.max([d for _, d in graph1.degree()])
  max_degree_g2 = np.max([d for _, d in graph2.degree()])

  # Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
  C = np.zeros((graph1.number_of_nodes(), graph2.number_of_nodes()))

  for i in range(graph1.number_of_nodes()):
    for j in range(graph2.number_of_nodes()):
      # Calculo de que tan cerca estan de tener la misma cantidad de conexicones que la mayor pareja
      temp_sum = float(graph1.degree[get_node_name_by_pos(graph1, i)] + graph2.degree[get_node_name_by_pos(graph2, j)])
      if (temp_sum) == (max_degree_g1 + max_degree_g2):
        nodes_sim_struct = 0.0
      else:
        nodes_sim_struct = 1.0 - (temp_sum / float(max_degree_g1 + max_degree_g2))

      # Calculo de que tan similar es la conexion entre vecinos
      total_nodes_sim_struct = structural_similarity_node(graph1, graph2,
                                                          get_node_name_by_pos(graph1, i),
                                                          get_node_name_by_pos(graph2, j))

      # Similitud entre los noso
      sim_node = similarity_node(get_node_by_position(graph1, i), get_node_by_position(graph2, j))

      C[i][j] = nodes_sim_struct + total_nodes_sim_struct + sim_node

  return C


def calculate_Z(graph1, graph2):
  # Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
  Z = np.zeros((graph1.number_of_nodes(), graph2.number_of_nodes()))

  # Cargar matrix
  for i in range(graph1.number_of_nodes()):
    for j in range(graph2.number_of_nodes()):
      Z[i][j] = similarity_node(get_node_by_position(graph1, i), get_node_by_position(graph2, j))

  return Z


def get_promising_neighbors(graph1, graph2, g1_mapping, g2_mapping, node1, node2, min_similarity_value, Z):
  # Vecinos de los nodos vistos
  total_g1_neighbors = []
  total_g2_neighbors = []

  # Lista temporal de pares de vecinos a revisar
  test_pairs = []

  # Se optienen los vecinos que no han sido mapeados ya
  g1_cycle = np.setdiff1d(list(graph1.neighbors(node1)), g1_mapping)
  g2_cycle = np.setdiff1d(list(graph2.neighbors(node2)), g2_mapping)

  for i in g1_cycle:
    for j in g2_cycle:
      pos_i = get_node_index(graph1, i)
      pos_j = get_node_index(graph2, j)
      test_pairs.append([Z[pos_i][pos_j], pos_i, pos_j])

      total_g1_neighbors.append(i)
      total_g2_neighbors.append(j)

  # Lista final de parejas
  final_result = []

  for i in test_pairs:
    if i[0] <= min_similarity_value:
      final_result.append([i[1], i[2]])

  return total_g1_neighbors, total_g2_neighbors, final_result


# Se setea una fila y columna dada con un valor dado
def set_matrix_ij_val(matrix, i, j, val):
  matrix[i, :] = val
  matrix[:, j] = val
  return matrix


# Dado dos nodos se retorna la cantidad de nodos compartidos, mapeados
def ok_neighbors(graph1, graph2, node1, node2, equivalences_list):
  test_pairs = []

  for i in graph1.neighbors(node1):
    for j in graph2.neighbors(node2):
      test_pairs.append([i, j])

  return len(intersect2d(test_pairs, equivalences_list))


# Calculo de matriz de nodos compartidos
def compute_neighbors_matrix(graph1, graph2, Z, equivalences_list):
  number_of_nodes1 = graph1.number_of_nodes()
  number_of_nodes2 = graph2.number_of_nodes()

  N = np.zeros((number_of_nodes1, number_of_nodes2))

  for i in range(number_of_nodes1):
    for j in range(number_of_nodes2):
      if Z[i][j] != np.infty:
        N[i][j] = ok_neighbors(graph1, graph2,
                               get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j),
                               equivalences_list)

  return N


def get_best_ij(N, Z, max_similarity_value):
  i_list, j_list = np.where(N == max_similarity_value)

  # Se colocan los valores por default
  actual_i = i_list[0]
  actual_j = j_list[0]
  actual_val = Z[actual_i][actual_j]

  for i in range(1, len(i_list)):
    x = Z[i_list[i]][j_list[i]]

    if actual_val > x:
      actual_val = x
      actual_i = i_list[i]
      actual_j = j_list[i]

    if actual_val == x:
      flat = random.choice([True, False])
      if flat:
        actual_val = x
        actual_i = i_list[i]
        actual_j = j_list[i]

  return actual_i, actual_j


def get_pos_value(matrix, value):
  pos_list = np.column_stack(np.where(matrix == value))
  result = random.choice(pos_list)
  return result[0], result[1]


def graph_aligning2(graph1, graph2, min_permit_value, repair_external_values=True):
  C = calculate_C(graph1, graph2)
  Z = calculate_Z(graph1, graph2)

  print(C)
  print(Z)

  result = []
  error_value_result = 0

  # Se almacenan todos los nodos que fueron mapeados (su valor real)
  g1_mapping = []
  g2_mapping = []

  most_similarity_value = np.amin(C)
  while most_similarity_value <= min_permit_value:
    i, j = get_pos_value(C, most_similarity_value)

    C, Z, error_value_result = expansion_cycle(C, Z,
                                               error_value_result,
                                               g1_mapping, g2_mapping,
                                               graph1, graph2,
                                               i, j,
                                               C[i][j],
                                               min_permit_value, result)

    most_similarity_value = np.amin(C)

  # Retornar resultado si no se quiere que se haga el proceso de forsar match
  if not repair_external_values:
    return error_value_result + (min(graph1.number_of_nodes(), graph2.number_of_nodes()) - len(result)), result

  # Se computa la matrix N de vecinos y se obtiene el elemento con mas nodos compartidos
  N = compute_neighbors_matrix(graph1, graph2, Z, result)
  max_similarity_value = np.amax(N)

  while max_similarity_value > 0:
    # De todos los posibles valores a hacer match, hacer match con el mejor
    i, j = get_best_ij(N, Z, max_similarity_value)

    C, Z, error_value_result = expansion_cycle(C, Z,
                                               error_value_result,
                                               g1_mapping, g2_mapping,
                                               graph1, graph2,
                                               i, j,
                                               Z[i][j],
                                               min_permit_value, result)

    # Se computa la matrix N de vecinos y se obtiene el elemento con mas nodos compartidos
    N = compute_neighbors_matrix(graph1, graph2, Z, result)
    max_similarity_value = np.amax(N)

  return error_value_result + (min(graph1.number_of_nodes(), graph2.number_of_nodes()) - len(result)), result


def expansion_cycle(C, Z, error_value_result, g1_mapping, g2_mapping, graph1, graph2, i, j, match_root_vale,
                    min_permit_value, result):
  C, Z, error_value_result, queue_neighbors = get_match_set_values(C, Z,
                                                                   error_value_result,
                                                                   g1_mapping, g2_mapping,
                                                                   graph1, graph2,
                                                                   i, j,
                                                                   min_permit_value,
                                                                   match_root_vale,
                                                                   result)

  # Se orden de menor a mayor, para crecer por vecinos mas prometedores
  queue_neighbors = sorted(queue_neighbors, key=lambda pair: Z[pair[0]][pair[1]], reverse=False)

  # Crecemos en ancho para cada un de los vecinos
  while queue_neighbors:
    pos_ij = queue_neighbors.pop(0)
    i, j = pos_ij[0], pos_ij[1]
    if Z[i][j] == np.Inf:
      continue

    C, Z, error_value_result, new_queue_neighbors = get_match_set_values(C, Z,
                                                                         error_value_result,
                                                                         g1_mapping, g2_mapping,
                                                                         graph1, graph2,
                                                                         i, j,
                                                                         min_permit_value,
                                                                         Z[i][j],
                                                                         result)

    queue_neighbors += setdiff2d(new_queue_neighbors, queue_neighbors)

    # Se orden de menor a mayor, para crecer por vecinos mas prometedores
    queue_neighbors = sorted(queue_neighbors, key=lambda pair: Z[pair[0]][pair[1]], reverse=False)

  return C, Z, error_value_result


def get_match_set_values(C, Z, error_value_result, g1_mapping, g2_mapping, graph1, graph2, i, j, min_permit_value,
                         most_similarity_value, result):
  # Se agrega el nuevo match mas el error
  error_value_result += most_similarity_value
  result.append([get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j)])

  # Se guardan los valores mapeados
  g1_mapping.append(get_node_name_by_pos(graph1, i))
  g2_mapping.append(get_node_name_by_pos(graph2, j))

  # Se llenan con infinito las filas y columnas que no pueden ser utilizadas
  C = set_matrix_ij_val(C, i, j, np.Inf)
  Z = set_matrix_ij_val(Z, i, j, np.Inf)

  # Obtienen los vecinos de los nodos mapeados que son mayores a un umbral
  g1_neighbors, g2_neighbors, queue_neighbors = get_promising_neighbors(graph1,
                                                                        graph2,
                                                                        g1_mapping,
                                                                        g2_mapping,
                                                                        get_node_name_by_pos(graph1, i),
                                                                        get_node_name_by_pos(graph2, j),
                                                                        min_permit_value, Z)

  return C, Z, error_value_result, queue_neighbors
