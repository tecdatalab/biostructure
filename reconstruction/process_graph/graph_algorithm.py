from operator import itemgetter

import numpy as np
import random
from sklearn.metrics import mean_squared_error
from scipy.optimize import linear_sum_assignment

from general_utils.graph_utils import get_node_by_position, get_node_name_by_pos, get_node_index
from general_utils.list_utils import setdiff2d, intersect2d


class Matching:

  def __init__(self):
    self.dic_g1g2_possible = {}

  # Obtiene la cantidad de vecinos para todos los vecinos de un nodo dado
  def cal_degree_non(self, graph, node):
    total_neighbors = 0
    for n in graph.neighbors(node):
      total_neighbors += graph.degree[n]
    return total_neighbors

  # Que tan similares son los nodos
  def possible_valid_error_neighbors(self, graph1, graph2, node1, node2, min_sin_value):
    graph1_dicc = {}
    graph2_dicc = {}

    matrix = []
    cont_i = 0

    if (len(list(graph1.neighbors(node1))) == 0) or (len(list(graph2.neighbors(node2))) == 0):
      return 0.0

    for i in graph1.neighbors(node1):
      row = []
      cont_j = 0
      graph1_dicc[cont_i] = i
      for j in graph2.neighbors(node2):
        row.append(self.similarity_node_data(graph1.nodes[i], graph2.nodes[j]))
        graph2_dicc[cont_j] = j
        cont_j += 1
      cont_i += 1
      matrix.append(row)

    row_ind, col_ind = linear_sum_assignment(matrix)

    error_in_struc_match = 0
    for i in range(len(row_ind)):
      error_in_struc_match += self.possible_valid_error(graph1, graph2,
                                                   graph1_dicc[row_ind[i]],
                                                   graph2_dicc[col_ind[i]],
                                                   min_sin_value)

    return error_in_struc_match

  # Que tan similares son los nodos mediantre el error minimo cuadrado
  def similarity_node_data(self, node1, node2):
    mse = mean_squared_error(node1['zd_descriptors'], node2['zd_descriptors'])

    return mse  # Definir como 0 el valor de igualdad total, entre mas mayor menos se parecen

  def similarity_node_struct(self, graph1, graph2, node1, node2):
    gra1_can_n = float(len(list(graph1.neighbors(node1))))
    gra2_can_n = float(len(list(graph2.neighbors(node2))))

    min_value = min(gra1_can_n, gra2_can_n)
    max_value = max(gra1_can_n, gra2_can_n)

    if max_value == 0:
      result = 0.0
    else:
      result = 1.0 - (min_value / max_value)

    return result

  # Que tan similares son los vecinos de los nodos
  def possible_valid_error(self, graph1, graph2, node1, node2, min_sin_value):

    if self.dic_g1g2_possible.get(str(node1) + str(node2)) != None:
      return self.dic_g1g2_possible.get(str(node1) + str(node2))

    if (len(list(graph1.neighbors(node1))) == 0) or (len(list(graph2.neighbors(node2))) == 0):
      self.dic_g1g2_possible[str(node1) + str(node2)] = 0.0
      return 0.0

    matrix = []
    for i in graph1.neighbors(node1):
      row = []
      for j in graph2.neighbors(node2):
        row.append(self.similarity_node_data(graph1.nodes[i], graph2.nodes[j]))
      matrix.append(row)

    row_ind, col_ind = linear_sum_assignment(matrix)

    ok_match = 0
    for i in range(len(row_ind)):
      if matrix[row_ind[i]][col_ind[i]] <= min_sin_value:
        ok_match += 1

    if ok_match == 0:
      result = 1.0
    else:
      result = (float(ok_match)) / float(min(len(list(graph1.neighbors(node1))), len(list(graph2.neighbors(node2)))))
    result = 1.0 - result

    self.dic_g1g2_possible[str(node1) + str(node2)] = result
    return result

  def calculate_C(self, graph1, graph2, min_similarity_value, g1_mapping=[], g2_mapping=[], Z=None):
    nodes_g1_check = np.setdiff1d(graph1.nodes, g1_mapping).tolist()
    nodes_g2_check = np.setdiff1d(graph2.nodes, g2_mapping).tolist()

    if nodes_g1_check == [] or nodes_g2_check == []:
      C = np.full((graph1.number_of_nodes(), graph2.number_of_nodes()), np.infty)
      return C

    max_degree_g1 = np.max([graph1.degree[node] for node in nodes_g1_check])
    max_degree_g2 = np.max([graph2.degree[node] for node in nodes_g2_check])

    mm_global_degree = min(max_degree_g1, max_degree_g2)

    # Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
    C = np.zeros((graph1.number_of_nodes(), graph2.number_of_nodes()))

    for i in range(graph1.number_of_nodes()):
      for j in range(graph2.number_of_nodes()):
        if Z is not None:
          if Z[i][j] == np.infty:
            C[i][j] = np.infty
            continue
        # Calculo de que tan cerca estan de tener la misma cantidad de conexicones que la mayor pareja
        temp_sum = float(min(graph1.degree[get_node_name_by_pos(graph1, i)], mm_global_degree)
                         + min(graph2.degree[get_node_name_by_pos(graph2, j)], mm_global_degree))

        if (temp_sum) == (mm_global_degree * 2):
          nodes_sim_root_struct = 0.0
        else:
          nodes_sim_root_struct = 1.0 - (temp_sum / float(mm_global_degree * 2))

        # Calculo de que tan similar es la conexion entre vecinos
        total_nodes_sim_struct = self.possible_valid_error_neighbors(graph1, graph2,
                                                                     get_node_name_by_pos(graph1, i),
                                                                     get_node_name_by_pos(graph2, j),
                                                                     min_similarity_value)

        # Similitud entre los nodo
        sim_node_data = self.similarity_node_data(get_node_by_position(graph1, i), get_node_by_position(graph2, j))

        # Ajuste de vecino
        insertion_price = self.possible_valid_error(graph1, graph2,
                                                    get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j),
                                                    min_similarity_value)

        C[i][j] = nodes_sim_root_struct + sim_node_data + total_nodes_sim_struct + insertion_price

    return C

  def calculate_Z(self, graph1, graph2, min_similarity_value):
    # Crear matrix n(cantidad de aristas g1) x m(cantidad de aristas g2)
    Z = np.zeros((graph1.number_of_nodes(), graph2.number_of_nodes()))

    # Cargar matrix
    for i in range(graph1.number_of_nodes()):
      for j in range(graph2.number_of_nodes()):
        # Ajuste de vecino
        insertion_price = self.possible_valid_error(graph1, graph2,
                                                    get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j),
                                                    min_similarity_value)

        # Similitud entre los nodo
        sim_node = self.similarity_node_data(get_node_by_position(graph1, i), get_node_by_position(graph2, j))

        Z[i][j] = sim_node + insertion_price

    return Z

  def get_promising_neighbors(self, graph1, graph2, g1_mapping, g2_mapping, node1, node2, min_similarity_value, Z):
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
  def set_matrix_ij_val(self, matrix, i, j, val):
    matrix[i, :] = val
    matrix[:, j] = val
    return matrix

  # Dado dos nodos se retorna la cantidad de nodos compartidos, mapeados
  def ok_neighbors(self, graph1, graph2, node1, node2, equivalences_list):
    test_pairs = []

    for i in graph1.neighbors(node1):
      for j in graph2.neighbors(node2):
        test_pairs.append([i, j])

    return len(intersect2d(test_pairs, equivalences_list))

  # Calculo de matriz de nodos compartidos
  def compute_neighbors_matrix(self, graph1, graph2, Z, equivalences_list):
    number_of_nodes1 = graph1.number_of_nodes()
    number_of_nodes2 = graph2.number_of_nodes()

    N = np.zeros((number_of_nodes1, number_of_nodes2))

    for i in range(number_of_nodes1):
      for j in range(number_of_nodes2):
        if Z[i][j] != np.infty:
          N[i][j] = self.ok_neighbors(graph1, graph2,
                                      get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j),
                                      equivalences_list)

    return N

  def get_best_ij(self, N, Z, max_similarity_value):
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

  def get_pos_value(self, matrix, value):
    pos_list = np.column_stack(np.where(matrix == value))
    result = random.choice(pos_list)
    return result[0], result[1]

  def do_graph_aligning(self, graph1, graph2, min_permit_value, repair_external_values=True):
    C = self.calculate_C(graph1, graph2, min_permit_value)
    Z = self.calculate_Z(graph1, graph2, min_permit_value)

    result = []
    error_value_result = 0

    # Se almacenan todos los nodos que fueron mapeados (su valor real)
    g1_mapping = []
    g2_mapping = []

    most_similarity_value = np.amin(C)
    while most_similarity_value <= min_permit_value:
      i, j = self.get_pos_value(C, most_similarity_value)

      C, Z, error_value_result = self.expansion_cycle(C, Z,
                                                 error_value_result,
                                                 g1_mapping, g2_mapping,
                                                 graph1, graph2,
                                                 i, j,
                                                 C[i][j],
                                                 min_permit_value, result)

      C = self.calculate_C(graph1, graph2, min_permit_value, g1_mapping, g2_mapping, Z)
      most_similarity_value = np.amin(C)

    # Retornar resultado si no se quiere que se haga el proceso de forsar match
    if not repair_external_values:
      if len(result) == 0:
        return error_value_result, result
      else:
        return error_value_result / len(result), result

    # Se computa la matrix N de vecinos y se obtiene el elemento con mas nodos compartidos
    N = self.compute_neighbors_matrix(graph1, graph2, Z, result)
    max_similarity_value = np.amax(N)

    while max_similarity_value > 0:
      # De todos los posibles valores a hacer match, hacer match con el mejor
      i, j = self.get_best_ij(N, Z, max_similarity_value)

      C, Z, error_value_result = self.expansion_cycle(C, Z,
                                                 error_value_result,
                                                 g1_mapping, g2_mapping,
                                                 graph1, graph2,
                                                 i, j,
                                                 Z[i][j],
                                                 min_permit_value, result)

      # Se computa la matrix N de vecinos y se obtiene el elemento con mas nodos compartidos
      N = self.compute_neighbors_matrix(graph1, graph2, Z, result)
      max_similarity_value = np.amax(N)

    if len(result) == 0:
      return error_value_result, result
    else:
      return error_value_result / len(result), result

  def expansion_cycle(self, C, Z, error_value_result, g1_mapping, g2_mapping, graph1, graph2, i, j, match_root_vale,
                      min_permit_value, result):
    C, Z, error_value_result, queue_neighbors = self.get_match_set_values(C, Z,
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

      C, Z, error_value_result, new_queue_neighbors = self.get_match_set_values(C, Z,
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


  def get_match_set_values(self, C, Z, error_value_result, g1_mapping, g2_mapping, graph1, graph2, i, j, min_permit_value,
                           most_similarity_value, result):
    # Se calcula el error estructural que conyeba este match
    error_node_struct = self.similarity_node_struct(graph1, graph2,
                                               get_node_name_by_pos(graph1, i),
                                               get_node_name_by_pos(graph2, j))

    # Se agrega el nuevo match mas el error
    error_value_result += most_similarity_value + error_node_struct
    result.append([get_node_name_by_pos(graph1, i), get_node_name_by_pos(graph2, j)])

    # Se guardan los valores mapeados
    g1_mapping.append(get_node_name_by_pos(graph1, i))
    g2_mapping.append(get_node_name_by_pos(graph2, j))

    # Se llenan con infinito las filas y columnas que no pueden ser utilizadas
    C = self.set_matrix_ij_val(C, i, j, np.Inf)
    Z = self.set_matrix_ij_val(Z, i, j, np.Inf)

    # Obtienen los vecinos de los nodos mapeados que son mayores a un umbral
    g1_neighbors, g2_neighbors, queue_neighbors = self.get_promising_neighbors(graph1,
                                                                          graph2,
                                                                          g1_mapping,
                                                                          g2_mapping,
                                                                          get_node_name_by_pos(graph1, i),
                                                                          get_node_name_by_pos(graph2, j),
                                                                          min_permit_value, Z)

    return C, Z, error_value_result, queue_neighbors


def graph_aligning(graph1, graph2, min_permit_value, repair_external_values=True):
  algoritm = Matching()
  return algoritm.do_graph_aligning(graph1, graph2, min_permit_value, repair_external_values)
