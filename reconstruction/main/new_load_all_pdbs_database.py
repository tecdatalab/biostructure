import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from general_utils.database_utils import get_chains_pdb_db
import numpy as np
import os
import random
import general_utils
import pandas as pd
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from general_utils.workspace_utils import is_work_in_cluster

from general_utils.pdb_utils import get_all_pdb_work, get_pdb_adn_arn, get_percentage_pbs_check_file
from general_utils.database_utils import get_chains_pdb_db, get_all_archive_pdb
from general_utils.temp_utils import clean_work_dir


def add_pdb(pdb_name):
  print("Start adding", pdb_name, flush=True)
  get_chains_pdb_db(pdb_name)

  print("Added", pdb_name, flush=True)
  with open("./all_check.txt", 'a+') as f:
    f.write(str(pdb_name) + "\n")


def get_to_load():
  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'
  know_pdb_path_cif_flag = os.path.dirname(__file__) + '/../files/pdb_can_chains_with_pdb_flag.csv'

  # All pdbs
  to_load = []
  if os.path.exists(know_pdb_path):
    pd_data_frame = pd.read_csv(know_pdb_path)
    for i in pd_data_frame.values.tolist():
      if i[1] == 1:
        to_load.append(i[0])

  # All pdbs not se cif file
  to_load_ok = []
  if os.path.exists(know_pdb_path_cif_flag):
    pd_data_frame = pd.read_csv(know_pdb_path_cif_flag)
    for i in pd_data_frame.values.tolist():
      if i[3] == True:
        to_load_ok.append(i[0])

  to_load = np.setdiff1d(np.array(to_load), np.array(get_all_archive_pdb())).tolist()
  to_load = np.setdiff1d(np.array(to_load), np.array(get_pdb_adn_arn())).tolist()
  to_load = np.intersect1d(np.array(to_load), np.array(to_load_ok)).tolist()

  return to_load


def gen_load_database():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      #all_names = get_to_load()
      all_names = ["1yfq"]

      random.shuffle(all_names)
      print("To load", len(all_names), flush=True)

      parallel_jobs = []
      for pdb_name in all_names:
        parallel_jobs.append([pdb_name, executor.submit(add_pdb, pdb_name)])

      total_to_do = len(all_names)
      actual_do = 0
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          print(f[0], e, flush=True)

        actual_do += 1
        print("Done", total_to_do, actual_do, actual_do / total_to_do, flush=True)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_load_database_file"
  else:
    general_utils.temp_utils.global_temp_dir = None
  gen_load_database()
