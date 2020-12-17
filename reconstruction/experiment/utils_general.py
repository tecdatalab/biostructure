import os
import pickle
import shutil
import tempfile
import random

from general_utils.download_utils import get_all_pdb_name, download_pdb
from pandas.core.common import flatten
from math import ceil
from pdb_to_mrc.miscellaneous import get_chains


def pdb_percentage(percentage):
  # Add ignore files
  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_no_work.txt'
  f_evil_pdb = open(evil_pdb_path)
  ignore_pdbs = f_evil_pdb.read().splitlines()

  # List know can chains pdb
  know_can_chains_pdb_path = os.path.dirname(__file__) + '/../files/pdb_can_chains.pkl'

  if os.path.exists(know_can_chains_pdb_path):
    can_chains_pdb_file = open(know_can_chains_pdb_path, "rb")
    can_chains_list = pickle.load(can_chains_pdb_file)
    can_chains_pdb_file.close()
    can_chains_list_var = list(flatten(can_chains_list))
  else:
    can_chains_list_var = []
    can_chains_list = []

  # Process
  all_names = get_all_pdb_name()  # 169315
  all_names = [all_names[0], all_names[1], all_names[3]]
  for i in range(len(all_names) - 1, -1, -1):
    if all_names[i] in ignore_pdbs or all_names[i] in can_chains_list_var:
      all_names.pop(i)

  # Add chains
  dirpath = tempfile.mkdtemp()

  for pdb_name in all_names:
    download_pdb(pdb_name, '{0}/{1}.pdb'.format(dirpath, pdb_name))
    chains = get_chains('{0}/{1}.pdb'.format(dirpath, pdb_name))
    can_chains_list.append([pdb_name, len(chains)])
    os.remove('{0}/{1}.pdb'.format(dirpath, pdb_name))

  shutil.rmtree(dirpath)
  if os.path.exists(know_can_chains_pdb_path):
    os.remove(know_can_chains_pdb_path)

  can_chains_pdb_file = open(know_can_chains_pdb_path, "wb")
  pickle.dump(can_chains_list, can_chains_pdb_file)
  can_chains_pdb_file.close()

  # Dic with can chains to pdbs
  dic_chains = {}
  for i in can_chains_list:
    can = i[1]
    pdb_name = i[0]

    if can in dic_chains:
      dic_chains[can].append(pdb_name)

    else:
      dic_chains[can] = [pdb_name]

  # Gen result
  result = []
  percentage_value = percentage / 100

  for i in dic_chains.keys():
    real_can_add = ceil(len(dic_chains[i]) * percentage_value)
    random.shuffle(dic_chains[i])
    result += dic_chains[i][:real_can_add]

  return result


def remove_get_dirs(path, force=False):
  result = []
  complete_path = os.path.abspath(path)
  list_dirs = os.listdir(complete_path)

  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_no_work.txt'
  f_evil_pdb = open(evil_pdb_path, 'a+')

  for dir_name in list_dirs:
    check_path = '{0}/{1}'.format(complete_path, dir_name)
    files_dir = os.listdir(check_path)

    if len(files_dir) == 1 and files_dir[0].split('.')[1] == 'csv':
      result.append(dir_name)
    else:

      if force:
        all_pdb = True
        for i in files_dir:
          if i.split('.')[1] != 'pdb':
            all_pdb = False

        if all_pdb:
          f_evil_pdb.write(dir_name + '\n')

      shutil.rmtree(check_path)

  f_evil_pdb.close()
  return result
