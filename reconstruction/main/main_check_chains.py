import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import numpy as np
import general_utils
from general_utils.cif_utils import get_chains_cif, cif_to_pdb
from general_utils.database_utils import get_all_archive_pdb, get_chains_pdb_db, delete_pdb_db
from general_utils.download_utils import download_pdb, download_cif
import os

from general_utils.temp_utils import gen_dir, free_dir
from general_utils.workspace_utils import is_work_in_cluster


def get_chains(pdb_name):
  dir_path = gen_dir()

  from general_utils.pdb_utils import get_chains_pdb

  flag_PDB = True
  try:
    path_of_pdb = '{0}/{1}.pdb'.format(dir_path, pdb_name)
    download_pdb(pdb_name, path_of_pdb)
    chains = get_chains_pdb(path_of_pdb)
  except:
    flag_PDB = False
    path_of_cif = '{0}/{1}.cif'.format(dir_path, pdb_name)
    download_cif(pdb_name, path_of_cif)
    chains = get_chains_cif(path_of_cif)

  free_dir(dir_path)

  return flag_PDB, chains


def check_pdb(pdb_name):
  # print("Enter", pdb_name, flush=True)
  flag_PDB, chains = get_chains(pdb_name)
  if (not flag_PDB) or (chains != get_chains_pdb_db(pdb_name)):
    delete_pdb_db(pdb_name)
    print("To update", pdb_name, flush=True)
    # get_chains_pdb_db(pdb_name)
    # print("Update", pdb_name, flush=True)
  else:
    with open("./all_check.txt", 'a+') as f:
      f.write(str(pdb_name) +"\n")


def check_chains():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      all_names = get_all_archive_pdb()
      total_to_do = len(all_names)

      if os.path.exists("./all_check.txt"):
        with open("./all_check.txt") as f:
          lines = f.read().splitlines()
          all_names = np.setdiff1d(np.array(all_names), np.array(lines)).tolist()

      total_to_do = len(all_names)
      parallel_jobs = []
      actual_do = 0
      for pdb_name in all_names:
        parallel_jobs.append([pdb_name, executor.submit(check_pdb, pdb_name)])
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          print(f[0], e, flush=True)

        actual_do += 1
        print("Done", f[0], total_to_do, actual_do, actual_do / total_to_do, flush=True)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_check_chains"
  else:
    general_utils.temp_utils.global_temp_dir = None
  # clean_work_dir()
  check_chains()
