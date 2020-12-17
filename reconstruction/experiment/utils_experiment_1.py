from math import ceil
import threading

from sklearn.neighbors import KDTree

from process_graph.process_segment_faces import get_n_points_cube


def gen_thread_points(id, mask, n_points_face, filter_value, semaphore, dic_result):
  result = get_n_points_cube(mask, n_points_face, filter_value)
  semaphore.acquire()
  dic_result[id] = result
  semaphore.release()


def gen_thread_KDTree(id, data, semaphore, dic_result):
  result = KDTree(data)
  semaphore.acquire()
  dic_result[id] = result
  semaphore.release()


def gen_keys_experiemnts(segments, n_points_face, filter_value, point_test, percentage_segments):
  face_points = {}
  kd_trees = {}

  sem = threading.Semaphore()
  thread_list = []
  for i in segments:
    t = threading.Thread(target=gen_thread_points, args=[i.id_segment, i.mask, n_points_face, filter_value,
                                                         sem, face_points])

    thread_list.append(t)
    t.start()

  for t in thread_list:
    t.join()

  thread_list = []
  for i in segments:
    t = threading.Thread(target=gen_thread_KDTree, args=[i.id_segment, face_points[i.id_segment], sem,
                                                         kd_trees])
    thread_list.append(t)
    t.start()

  for t in thread_list:
    t.join()

  result = []
  for segment in segments:
    # check code
    dist, _ind = kd_trees[segment.id_segment].query([point_test], k=1)
    result.append([segment.id_segment, dist[0][0]])

  # print(result)
  result = sorted(result, key=lambda val: val[1])
  # print(result)

  final_result = [i[0] for i in result[:int(ceil(len(result) * ((100 - percentage_segments) / 100)))]]
  # print(final_result)
  final_result.sort()
  # print(final_result)
  return final_result
