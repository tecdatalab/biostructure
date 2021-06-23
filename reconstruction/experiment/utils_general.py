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


def pdb_percentage(percentage, executor=None, min_can_chains=0):
  from general_utils.database_utils import get_chains_pdb_db

  # Process
  can_chains_list_dic = {}
  all_names = get_all_pdb_work()  # 169315
  for pdb_name in all_names:
    chains_len = len(get_chains_pdb_db(pdb_name))
    exist_flag = can_chains_list_dic.get(chains_len, None)
    if exist_flag is None:
      can_chains_list_dic[chains_len] = [pdb_name]
    else:
      can_chains_list_dic[chains_len].append(pdb_name)

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
