import os
import shutil
import random
from ast import literal_eval
import numpy as np
import pandas as pd
from csv_modules.csv_writer import write_in_file
from general_utils.cif_utils import get_chains_cif, cif_to_pdb
from general_utils.download_utils import download_pdb, download_cif
from math import ceil

from general_utils.pdb_utils import get_chains_pdb, get_cube_pdb, move_pdb_center, \
  get_all_pdb_work
from general_utils.temp_utils import gen_dir, free_dir


def get_parallel_can_chains_chunk(pdb_name, dir_path):
  from general_utils.pdb_utils import get_chains_pdb

  try:
    path_of_pdb = '{0}/{1}.pdb'.format(dir_path, pdb_name)
    download_pdb(pdb_name, path_of_pdb)
    chains = get_chains_pdb(path_of_pdb)
  except:
    path_of_cif = '{0}/{1}.cif'.format(dir_path, pdb_name)
    download_cif(pdb_name, path_of_cif)
    chains = get_chains_cif(path_of_cif)
    path_of_pdb = os.path.join(dir_path, "pdbFile.pdb")
    cif_to_pdb(path_of_cif, path_of_pdb)
    os.remove(path_of_cif)

  move_pdb_center(path_of_pdb)
  cube_dimensions = get_cube_pdb(path_of_pdb)

  result = [pdb_name, len(chains), cube_dimensions]
  os.remove(path_of_pdb)
  return result


def pdb_percentage(percentage, executor=None, min_can_chains=0):
  headers_csv = ['Pdb', 'Can Chains', 'Cube dimension']

  # List know can chains pdb
  know_can_chains_pdb_path = os.path.dirname(__file__) + '/../files/pdb_can_chains.csv'

  if os.path.exists(know_can_chains_pdb_path):
    pd_data_frame = pd.read_csv(know_can_chains_pdb_path,
                                converters={"Cube dimension": literal_eval})
    can_chains_list_name = pd_data_frame["Pdb"].tolist()
  else:
    can_chains_list_name = []

  # Process
  all_names = get_all_pdb_work()  # 169315
  # all_names = all_names[:3]
  all_names = np.setdiff1d(np.array(all_names), np.array(can_chains_list_name))
  all_names = all_names.tolist()

  # Add chains
  dirpath = gen_dir()

  flag_error = False

  if executor is not None:
    parallel_jobs = []
    for i in all_names:
      parallel_jobs.append([executor.submit(get_parallel_can_chains_chunk, i, dirpath)])

    for f in parallel_jobs:
      try:
        list_append = f[0].result()
        data_write = [list_append]
        write_in_file(know_can_chains_pdb_path, headers_csv, data_write)
      except Exception as e:
        flag_error = True
  else:
    for i in all_names:
      list_append = get_parallel_can_chains_chunk(i, dirpath)
      data_write = [list_append]
      write_in_file(know_can_chains_pdb_path, headers_csv, data_write)

  free_dir(dirpath)

  if flag_error:
    raise NameError('Error in download pdbs')

  # Gen data use
  pd_data_frame = pd.read_csv(know_can_chains_pdb_path,
                              converters={"Cube dimension": literal_eval})
  can_chains_list_dic = {}

  for value in pd_data_frame.values.tolist():
    exist_flag = can_chains_list_dic.get(value[1], None)
    if exist_flag == None:
      can_chains_list_dic[value[1]] = [value[0]]
    else:
      can_chains_list_dic[value[1]].append(value[0])

  # Gen result
  result = []
  percentage_value = percentage / 100

  for i in can_chains_list_dic.keys():
    if i >= min_can_chains:
      real_can_add = ceil(len(can_chains_list_dic[i]) * percentage_value)
      random.shuffle(can_chains_list_dic[i])
      result += can_chains_list_dic[i][:real_can_add]

  return result


def remove_get_dirs(path, add_to_ignore_files=False, can_csv=1):
  result = []
  complete_path = os.path.abspath(path)
  list_dirs = os.listdir(complete_path)

  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_no_work.txt'
  f_evil_pdb = open(evil_pdb_path, 'a+')

  for dir_name in list_dirs:
    check_path = '{0}/{1}'.format(complete_path, dir_name)
    files_dir = os.listdir(check_path)

    flag = True
    for i in files_dir:
      if i.find(".") == -1:
        flag = False
        break
      if not i.split('.')[1] == 'csv':
        flag = False
        break

    if len(files_dir) == can_csv and flag:
      result.append(dir_name)

    else:
      if add_to_ignore_files:
        all_pdb = True
        for i in files_dir:
          if i.find(".") == -1:
            all_pdb = False
            break
          if i.split('.')[1] != 'pdb':
            all_pdb = False
            break

        if all_pdb:
          f_evil_pdb.write(dir_name + '\n')

      shutil.rmtree(check_path)

  f_evil_pdb.close()
  return result
