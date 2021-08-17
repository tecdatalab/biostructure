import os
import shutil
import random
from ast import literal_eval
import numpy as np
import pandas as pd
import errno

from csv_modules.csv_writer import write_in_file
from general_utils.cif_utils import get_chains_cif, cif_to_pdb, cif_to_chains_pdb_files
from general_utils.download_utils import download_pdb, download_cif
from math import ceil
from pymol import cmd

from general_utils.pdb_utils import get_chains_pdb, get_cube_pdb, move_pdb_center, \
  get_all_pdb_work
from general_utils.temp_utils import gen_dir, free_dir


def pdb_percentage(percentage, executor=None, min_can_chains=0):
  from general_utils.database_utils import get_chains_pdb_db
  from general_utils.database_utils import get_dicc_pdbs_can_chains

  # Process
  can_chains_list_dic = get_dicc_pdbs_can_chains()

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


def make_dir_pdb(work_dir, pdb_id):
  complete_path = os.path.abspath(work_dir)
  dirs = os.listdir(complete_path)

  if pdb_id in dirs:
    if len(os.listdir(os.path.join(complete_path, pdb_id))) > 0:
      return
    else:
      os.rmdir(os.path.join(complete_path, pdb_id))

  from general_utils.database_utils import get_db_chains_files_db
  get_db_chains_files_db(pdb_id, work_dir)


def check_RMSD_result_algorithm(work_dir, all_chains, check_list, original_pdb, changed_pdb, changed_chain,
                                pdb_work_chain):
  final_result = {}
  final_result_list = []

  if not os.path.exists(work_dir):
    os.mkdir(work_dir)

  make_dir_pdb(work_dir, original_pdb)
  make_dir_pdb(work_dir, changed_pdb)

  for i in check_list:
    #Restar
    cmd.reinitialize()

    pdb_chaini = original_pdb
    pdb_chainj = original_pdb

    ipos = i[0]
    jpos = i[1]

    chain1 = all_chains[ipos - 1]
    chain2 = all_chains[jpos - 1]

    if chain2 == changed_chain:
      chain2 = pdb_work_chain
      pdb_chainj = changed_pdb

    path_chain_i = work_dir + "/" + pdb_chaini + "/" + "{}_{}.pdb".format(pdb_chaini, chain1)
    path_chain_j = work_dir + "/" + pdb_chainj + "/" + "{}_{}.pdb".format(pdb_chainj, chain2)

    namei = "{}_{}".format(pdb_chaini + "i", chain1)
    namej = "{}_{}".format(pdb_chainj + "j", chain2)

    if (not os.path.exists(path_chain_i)):
      # from general_utils.database_utils import delete_pdb_db
      # delete_pdb_db(pdb_chaini)
      # from general_utils.database_utils import get_chains_pdb_db
      # get_chains_pdb_db(pdb_chaini)
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path_chain_i)

    if (not os.path.exists(path_chain_j)):
      # from general_utils.database_utils import delete_pdb_db
      # delete_pdb_db(pdb_chainj)
      # from general_utils.database_utils import get_chains_pdb_db
      # get_chains_pdb_db(pdb_chainj)
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path_chain_j)

    cmd.load(path_chain_i, namei)
    cmd.load(path_chain_j, namej)

    try:
      result = cmd.align(namei, namej, cycles=1000)
      add_data = [cmd.count_atoms(namei), cmd.count_atoms(namej), result[0], result[1], result[5]]

    except Exception as e:
      print(f"Error in alinament {namei}, {namej}", flush=True)
      result = [9999999, 0, 0, 9999999, 0, 0, 0]
      add_data = [0, 0, result[0], result[1], result[5]]

      ##raise Exception(f"Error in alinament {namei}, {namej}").with_traceback(e.__traceback__)


    final_result[namei + "_" + namej] = add_data
    final_result_list.append(result[0])

  final_result["avg"] = averageOfList(final_result_list)
  return final_result


def averageOfList(num):
  sumOfNumbers = 0
  for t in num:
    sumOfNumbers = sumOfNumbers + t

  if len(num) == 0:
    return -1
  avg = sumOfNumbers / len(num)
  return avg


def check_RMSD_result_all(work_dir,
                          all_chains_A,
                          all_chains_B,
                          check_list,
                          pdb_A,
                          pdb_B):
  final_result = {}
  final_result_list = []

  if not os.path.exists(work_dir):
    os.mkdir(work_dir)

  make_dir_pdb(work_dir, pdb_A)
  make_dir_pdb(work_dir, pdb_B)

  for i in check_list:
    #Restar
    cmd.reinitialize()

    pdb_chaini = pdb_A
    pdb_chainj = pdb_B

    ipos = i[0]
    jpos = i[1]

    chain1 = all_chains_A[ipos - 1]
    chain2 = all_chains_B[jpos - 1]

    path_chain_i = work_dir + "/" + pdb_chaini + "/" + "{}_{}.pdb".format(pdb_chaini, chain1)
    path_chain_j = work_dir + "/" + pdb_chainj + "/" + "{}_{}.pdb".format(pdb_chainj, chain2)

    namei = "{}_{}".format(pdb_chaini + "i", chain1)
    namej = "{}_{}".format(pdb_chainj + "j", chain2)

    if (not os.path.exists(path_chain_i)):
      # from general_utils.database_utils import delete_pdb_db
      # delete_pdb_db(pdb_chaini)
      # from general_utils.database_utils import get_chains_pdb_db
      # get_chains_pdb_db(pdb_chaini)
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path_chain_i)

    if (not os.path.exists(path_chain_j)):
      # from general_utils.database_utils import delete_pdb_db
      # delete_pdb_db(pdb_chainj)
      # from general_utils.database_utils import get_chains_pdb_db
      # get_chains_pdb_db(pdb_chainj)
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path_chain_j)

    cmd.load(path_chain_i, namei)
    cmd.load(path_chain_j, namej)

    try:
      result = cmd.align(namei, namej)
      add_data = [cmd.count_atoms(namei), cmd.count_atoms(namej), result[0], result[1], result[5]]

    except Exception as e:
      print(f"Error in alinament {namei}, {namej}", flush=True)
      result = [9999999, 0, 0, 9999999, 0, 0, 0]
      add_data = [0, 0, result[0], result[1], result[5]]

      ##raise Exception(f"Error in alinament {namei}, {namej}").with_traceback(e.__traceback__)


    final_result[namei + "_" + namej] = add_data
    final_result_list.append(result[0])

  final_result["avg"] = averageOfList(final_result_list)
  return final_result
