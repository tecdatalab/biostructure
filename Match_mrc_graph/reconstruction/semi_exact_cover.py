import numpy as np
import copy
import threading

semaphore_run = None
semaphore_append = None
actual_thread_number = 0


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


def get_semi_exact_s_aux(matrix_p, actual_pos, max_len, result_list, max_thread_number):
    global semaphore_run, actual_thread_number, semaphore_append

    if actual_pos == max_len:
        semaphore_append.acquire()
        result_list.append([[count_set_bits(np.sum(matrix_p))] + matrix_p.tolist()])
        semaphore_append.release()
        return

    mask = get_binary_mask(actual_pos, 4)
    pos = np.where(np.bitwise_and(matrix_p, mask) == mask)[0]
    can_pos = pos.shape[0]
    me_threads = list()

    if can_pos != 0:
        for i in range(can_pos):
            delete_items = np.where(np.bitwise_and(matrix_p, matrix_p[pos[i]]) != 0)[0]
            delete_items = delete_items[delete_items != pos[i]]
            matrix_aux = np.delete(matrix_p, delete_items)
            # Thread
            get_semi_exact_s_aux_thread(actual_pos, matrix_aux, max_len, max_thread_number, me_threads, result_list)

        matrix_aux = np.delete(matrix_p, pos)
        if matrix_aux.size != 0:
            # Thread
            get_semi_exact_s_aux_thread(actual_pos, matrix_aux, max_len, max_thread_number, me_threads, result_list)
    else:
        # Thread
        get_semi_exact_s_aux_thread(actual_pos, matrix_p, max_len, max_thread_number, me_threads, result_list)

    for i in me_threads:
        i.join()
        semaphore_run.acquire()
        actual_thread_number -= 1
        semaphore_run.release()


def get_semi_exact_s_aux_thread(actual_pos, matrix_aux, max_len, max_thread_number, me_threads, result_list):
    global actual_thread_number, semaphore_run

    semaphore_run.acquire()
    if actual_thread_number < max_thread_number:
        thread = threading.Thread(target=get_semi_exact_s_aux, args=(matrix_aux, actual_pos + 1, max_len,
                                                                     result_list, max_thread_number))
        actual_thread_number += 1
        me_threads.append(thread)
        thread.start()
        semaphore_run.release()
    else:
        semaphore_run.release()
        get_semi_exact_s_aux(matrix_aux, actual_pos + 1, max_len, result_list, max_thread_number)


def get_semi_exact_s(matrix_p, max_len, top, max_threads=12):
    global semaphore_run, actual_thread_number, semaphore_append
    semaphore_run = threading.Semaphore()
    semaphore_append = threading.Semaphore()

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
    get_semi_exact_s_aux(matrix, 0, max_len, result, max_threads)
    # print(actual_thread_number)
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
