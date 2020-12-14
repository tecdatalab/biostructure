import os
import shutil
import itertools
import numpy as np
from sklearn.neighbors import KDTree
from process_graph.process_segment_faces import get_n_points_cube
import threading
import json
import requests


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

  final_result = [i[0] for i in result[:int(round(len(result) * ((100 - percentage_segments) / 100)))]]
  # print(final_result)
  final_result.sort()
  # print(final_result)
  return final_result


def remove_get_dirs(path):
  result = []
  complete_path = os.path.abspath(path)
  list_dirs = os.listdir(complete_path)

  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_no_work.txt'
  f_evil_pdb = open(evil_pdb_path, 'a+')

  list_dirs = sorted(list_dirs, key=lambda f: os.path.getmtime('{0}/{1}'.format(complete_path, f)))
  if list_dirs == []:
    return []
  max_modification_time = os.path.getmtime('{0}/{1}'.format(complete_path, list_dirs[-1]))

  for dir_name in list_dirs:
    check_path = '{0}/{1}'.format(complete_path, dir_name)
    actual_modification_time = os.path.getmtime(check_path)
    files_dir = os.listdir(check_path)

    if len(files_dir) == 1 and files_dir[0].split('.')[1] == 'csv':
      result.append(dir_name)
    else:

      all_pdb = True
      for i in files_dir:
        if i.split('.')[1] != 'pdb':
          all_pdb = False

      if all_pdb:
        f_evil_pdb.write(dir_name + '\n')

      shutil.rmtree(check_path)

  f_evil_pdb.close()
  return result


def get_similar_pdb(pdb_name):
  search_request = {
                    "query": {
                      "type": "terminal",
                      "service": "structure",
                      "parameters": {
                        "value": {
                          "entry_id": pdb_name.upper(),
                          "assembly_id": "1"
                        },
                        "operator": "relaxed_shape_match"
                      }
                    },
                    "return_type": "assembly",
                    "request_options": {
                      "pager": {
                        "start": 0,
                        "rows": 100
                      },
                      "scoring_strategy": "combined",
                      "sort": [
                        {
                          "sort_by": "score",
                          "direction": "desc"
                        }
                      ]
                    }
                  }
  json_dump = json.dumps(search_request)
  url_get = 'https://search.rcsb.org/rcsbsearch/v1/query?json={0}'.format(json_dump)
  response = requests.get(url_get)
  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"]:
    result.append([i["identifier"].split("-")[0].lower(), i["score"]])
  return result
