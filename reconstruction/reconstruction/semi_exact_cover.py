import copy

import numpy as np
import itertools
import threading

semaphore_run = None
semaphore_append = None
condition = None
current_workers = 0


def get_key_value(data):
  return data[0]


def create_result(list_key, dic, result_global, cover, actual_result=[], count_update=[]):
  if len(list_key) == 0:
    result_global.append([cover, actual_result])
    ### Count updates
    if count_update != []:
      count_update[0] += 1
    return

  for i in dic[list_key[0]]:
    create_result(list_key[1:], dic, result_global, cover, actual_result + [i], count_update=count_update)


# Function determines the neighbors of a given vertex
def neighbors(vertex, graph):
  return np.where(graph[vertex] == 1)[0]


# The Bron-Kerbosch algorithm
def bron_kerbosch_thread(work_list, result, graph, matrix_p, count_update=[]):
  global semaphore_run, semaphore_append, condition, current_workers
  while True:
    semaphore_run.acquire()
    if len(work_list) > 0:
      r, p, x = work_list.pop(0)
      current_workers -= 1
      ### Count updates
      if count_update != []:
        count_update[0] += 1
      # print(threading.get_ident(), "work")
    else:
      if current_workers == 0 and len(work_list) == 0:
        semaphore_run.release()
        return

      condition.acquire()
      semaphore_run.release()
      condition.wait()
      condition.release()

      semaphore_run.acquire()
      if current_workers == 0 and len(work_list) == 0:
        semaphore_run.release()
        return
      else:
        semaphore_run.release()
        continue
    semaphore_run.release()

    if len(p) == 0 and len(x) == 0:
      semaphore_run.acquire()
      current_workers += 1
      semaphore_run.release()
      continue

    for vertex in p[:]:
      r_new = r[::]
      r_new.append(vertex)
      all_neighbors = neighbors(vertex, graph)
      p_new = np.intersect1d(p, all_neighbors)  # p intersects N(vertex)
      x_new = np.intersect1d(x, all_neighbors)  # x intersects N(vertex)
      can_ones = bin(np.sum(matrix_p[r_new])).count("1")

      # Add new result
      semaphore_append.acquire()
      result.append([can_ones, matrix_p[r_new]])
      semaphore_append.release()

      # Add new possibility
      semaphore_run.acquire()

      ### Count updates
      if count_update != []:
        count_update[0] += 5 + 2 + len(np.where(p == vertex)[0])

      work_list.append([r_new, p_new, x_new])
      semaphore_run.release()

      p = np.delete(p, np.where(p == vertex))
      x = np.append(x, [vertex])

    semaphore_run.acquire()
    current_workers += 1
    condition.acquire()
    condition.notifyAll()
    condition.release()
    semaphore_run.release()


# The Bron-Kerbosch algorithm
def bron_kerbosch(r, p, x, result, graph, matrix_p, max_thread_number, count_update=[]):
  work_list = [[r, p, x]]
  me_threads = list()

  if count_update != []:
    count_update[0] += 1

  for i in range(max_thread_number):
    thread = threading.Thread(target=bron_kerbosch_thread, args=(work_list, result, graph, matrix_p, count_update,))
    me_threads.append(thread)
    thread.start()

  for i in me_threads:
    i.join()


def get_semi_exact_s_aux(matrix_p, result_list, max_thread_number, count_update=[]):
  can_candidates = matrix_p.shape[0]
  adjacency_matrix = np.zeros((can_candidates, can_candidates), dtype=np.bool)
  ### Count updates
  if count_update != []:
    count_update[0] += 1

  for item in itertools.combinations(np.arange(can_candidates), 2):
    if np.bitwise_and(matrix_p[item[0]], matrix_p[[item[1]]]) == 0:
      adjacency_matrix[item[0]][item[1]] = 1
      adjacency_matrix[item[1]][item[0]] = 1
      ### Count updates
      if count_update != []:
        count_update[0] += 2

  bron_kerbosch([], np.arange(can_candidates), [], result_list, adjacency_matrix, matrix_p, max_thread_number,
                count_update=count_update)


def get_semi_exact_s(matrix, top, max_threads=12, count_update=[]):
  global semaphore_run, semaphore_append, condition
  semaphore_run = threading.Semaphore()
  semaphore_append = threading.Semaphore()
  condition = threading.Condition()
  matrix = copy.copy(matrix)
  ### Count updates
  if count_update != []:
    count_update[0] += 1

  dic_val_num = {}
  len_matrix = matrix.shape[0]

  for i in range(len_matrix - 1, -1, -1):
    if matrix[i] not in dic_val_num:
      dic_val_num[matrix[i]] = [i]
      ### Count updates
      if count_update != []:
        count_update[0] += 1
    else:
      dic_val_num[matrix[i]].append(i)
      matrix = np.delete(matrix, i)
      ### Count updates
      if count_update != []:
        count_update[0] += 2

  result = []
  get_semi_exact_s_aux(matrix, result, max_threads, count_update=count_update)
  result = sorted(result, key=get_key_value)  # Change

  ### Count updates
  if count_update != []:
    count_update[0] += len(result) * int(np.log(len(result)))

  combinations = []
  len_result = len(result)
  for i in range(len_result - 1, -1, -1):
    create_result(result[i][1], dic_val_num, combinations, result[i][0], count_update=count_update)

  return combinations[0:min(top, len(combinations))]
