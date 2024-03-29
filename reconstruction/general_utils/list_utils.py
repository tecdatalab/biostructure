import math
import random
import numpy as np
import itertools


def get_element_list(index, list_elements):
  result = []
  for i in list_elements:
    result.append(i[index])
  return result


def generate_binary_matrix(matrix):
  result = []
  for i in matrix:
    val = '0b' + ''.join([str(int(x)) for x in i])
    result.append(int(val, 2))
  result = np.array(result)

  return result


def combinations_12n(n):
  stuff = [i for i in range(n)]
  result = []

  for L in range(0, len(stuff) + 1):
    for subset in itertools.combinations(stuff, L):
      result.append(list(subset))
  return result


def combinations_i2jInK(n, min_can, max_can, max_combinations=1000, max_combinations_peer_value=500, check_max=True):
  stuff = [i for i in range(n)]
  result = []

  for k in range(min_can, max_can + 1):

    if not check_max or (combinations_formule(k,n) < max_combinations):
      for subset in itertools.combinations(stuff, k):
        result.append(list(subset))
    else:
      result += seudo_combinations_i2jInK(n, k, max_combinations_peer_value)
  return result

def combinations_formule(k, n):
  combi = math.comb(n, k)
  return combi


def seudo_combinations_i2jInK(n, can_add,  can_per_val):
  random.seed()

  values = [i for i in range(n)]
  result = []
  adding = 0

  while adding < can_per_val:
    temp_list = []
    random.shuffle(values)

    for actual_pos in range(n):
      if (can_add - len(temp_list)) == (n - actual_pos):
        temp_list += values[actual_pos:]

        temp_list.sort()
        if temp_list not in result:
          result.append(temp_list)
          adding += 1
        break

      if random.random() < .5:
        temp_list.append(values[actual_pos])

      if can_add == len(temp_list):
        temp_list.sort()
        if temp_list not in result:
          result.append(temp_list)
          adding += 1
        break

  return result


def intersect2d(listA, listB):
  a1_rows = set(map(tuple, listA))
  a2_rows = set(map(tuple, listB))
  result = [list(elem) for elem in list(a1_rows.intersection(a2_rows))]
  return result


def setdiff2d(list_a, list_b):
  a1_rows = set(map(tuple, list_a))
  a2_rows = set(map(tuple, list_b))
  result = [list(elem) for elem in list(a1_rows.difference(a2_rows))]
  return result
