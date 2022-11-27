import csv  # imports module csv
import os
import sys
import pathlib
import shutil
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

# from general_utils.database_utils import clear_collection

def remove_folder_parallel(folder_delete):
  shutil.rmtree(folder_delete, ignore_errors=True)


def main():

  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      folders_to_delete = [
        "/work/lcastillo/allCIFS_tmp",
        "/work/lcastillo/allPDBS_tmp",
        "/work/lcastillo/allPDBS_unzip",
        "/work/lcastillo/chain_pdb",
        "/work/lcastillo/allCIFS_unzip",
        "/work/lcastillo/all_pdb_files"
      ]

      #Create Works
      parallel_jobs = []
      for folder_to_delete in folders_to_delete:
          parallel_jobs.append(executor.submit(remove_folder_parallel, folder_to_delete))

      for f in parallel_jobs:
        try:
          result = f.result()
        except Exception as e:
          print(e, flush=True)

      print("Start clean database")
      # clear_collection()


main()
