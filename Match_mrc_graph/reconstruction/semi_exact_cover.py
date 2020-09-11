import numpy as np
import copy


def count_set_bits(n):
    count = 0
    while n:
        count += n & 1
        n >>= 1
    return count


def get_binary_mask(true_pos, can_elements):
    result_a = '0' * true_pos
    result_b = '0' * (can_elements - 1 - true_pos)
    result = result_a + '1' + result_b
    return int('0b' + result, 2)


def get_key_value(data):
    return data[0][0]


def create_result(list_key, dic, result_global, cover, actual_result=[]):
    if list_key == []:
        result_global.append([cover, actual_result])
        return

    for i in dic[list_key[0]]:
        create_result(list_key[1:], dic, result_global, cover, actual_result + [i])


def get_semi_exact_s_aux(matrix_p, actual_pos, max_len, result_list):
    if actual_pos == max_len:
        result_list.append([[count_set_bits(np.sum(matrix_p))] + matrix_p.tolist()])
        return

    mask = get_binary_mask(actual_pos, 4)
    pos = np.where(np.bitwise_and(matrix_p, mask) == mask)[0]
    can_pos = pos.shape[0]

    if can_pos != 0:
        for i in range(can_pos):
            delete_items = np.where(np.bitwise_and(matrix_p, matrix_p[pos[i]]) != 0)[0]
            delete_items = delete_items[delete_items != pos[i]]
            matrix_aux = np.delete(matrix_p, delete_items)
            get_semi_exact_s_aux(matrix_aux, actual_pos + 1, max_len, result_list)

        matrix_aux = np.delete(matrix_p, pos)
        if matrix_aux.size != 0:
            get_semi_exact_s_aux(matrix_aux, actual_pos + 1, max_len, result_list)
    else:
        get_semi_exact_s_aux(matrix_p, actual_pos + 1, max_len, result_list)


def get_semi_exact_s(matrix_p, max_len, top):
    matrix = copy.copy(matrix_p)
    dic_val_num = {}
    len_matrix = matrix.shape[0]

    for i in range(len_matrix - 1, -1, -1):
        if matrix[i] not in dic_val_num:
            dic_val_num[matrix[i]] = [i]
        else:
            dic_val_num[matrix[i]].append(i)
            matrix = np.delete(matrix, i)

    result = []
    get_semi_exact_s_aux(matrix, 0, max_len, result)
    result = np.array(result)
    result = np.unique(result)
    # print("Final result:", result)

    combinations = []
    len_result = len(result)
    for i in range(len_result - 1, -1, -1):
        create_result(result[i][1:], dic_val_num, combinations, result[i][0])

        if len(combinations) >= top:
            return combinations[0:top]

    return combinations
