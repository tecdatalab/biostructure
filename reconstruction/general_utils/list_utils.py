import numpy as np


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
  import itertools

  stuff = [i for i in range(n)]
  result = []

  for L in range(0, len(stuff) + 1):
    for subset in itertools.combinations(stuff, L):
      result.append(list(subset))
  return result


def union_multidimensional(listA, listB):
  A = np.array(listA)
  B = np.array(listB)

  nrows, ncols = A.shape
  dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
           'formats': ncols * [A.dtype]}

  C = np.intersect1d(A.view(dtype), B.view(dtype))

  # This last bit is optional if you're okay with "C" being a structured array...
  C = C.view(A.dtype).reshape(-1, ncols)
  return C


def diff_multidimensional(list_a, list_b):
  a1_rows = set(map(tuple, list_a))
  a2_rows = set(map(tuple, list_b))
  result = [list(elem) for elem in list(a1_rows.difference(a2_rows))]
  return result
