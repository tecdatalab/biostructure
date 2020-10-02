import numpy as np


def get_element_list(index, list_elements):
    result = []
    for i in list_elements:
        result.append(i[index])
    return result


def generate_binary_matrix(matrix):
    result = []
    for i in matrix:
        val = '0b' + ''.join([str(x) for x in i])
        result.append(int(val, 2))
    result = np.array(result)
    return result
