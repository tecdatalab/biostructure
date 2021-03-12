import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import numpy as np
import random
import general_utils
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from general_utils.workspace_utils import is_work_in_cluster

from general_utils.pdb_utils import get_all_pdb_work
from general_utils.database_utils import get_chains_pdb_db, get_all_archive_pdb
from general_utils.temp_utils import clean_work_dir


def add_pdb(pdb_name):
  get_chains_pdb_db(pdb_name)


def gen_load_database():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      all_names = np.setdiff1d(np.array(get_all_pdb_work()), np.array(get_all_archive_pdb())).tolist()
      random.shuffle(all_names)

      parallel_jobs = []
      for pdb_name in all_names:
        parallel_jobs.append([pdb_name, executor.submit(add_pdb, pdb_name)])
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          print(f[0])
          print(e)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_load_database"
  else:
    general_utils.temp_utils.global_temp_dir = None
  # clean_work_dir()
  gen_load_database()
