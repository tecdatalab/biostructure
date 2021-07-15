import sys
import pathlib


sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
from experiment.utils_general import check_RMSD_result_algorithm, make_dir_pdb
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import pandas as pd
from ast import literal_eval
import os


def checkFile(work_dir, pdb_name):
  pdb_name = pdb_name.replace("+", "")
  pdb_name = pdb_name.replace("5000000000000.0", "5e12")

  make_dir_pdb(work_dir, pdb_name)
  return True

def main():
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      file_name = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/data_experiment_1_b_v1_exe_1/4u0d/result_chain.csv"
      work_dir = "/work/lcastillo/RMSD"
      work_dir = "./RMSD"

      if not os.path.exists(work_dir):
        os.mkdir(work_dir)

      csv_file = pd.read_csv(file_name, converters={"Chains": literal_eval,
                                                    "Point Original": literal_eval,
                                                    "Point Test": literal_eval,
                                                    "Point Original syn": literal_eval,
                                                    "Point Test syn": literal_eval,
                                                    "Match": literal_eval
                                                    })

      col_names_dowloand = csv_file['Pdb'].tolist()
      col_names_dowloand += csv_file['Pdb work'].tolist()
      col_names_dowloand = np.unique(col_names_dowloand).tolist()

      parallel_dow = []
      for i in col_names_dowloand:
        parallel_dow.append([executor.submit(checkFile,
                                             work_dir,
                                             i),
                             i])

      for i in range(len(parallel_dow)):
        print(parallel_dow[i][1], i/len(parallel_dow), flush=True)
        print(parallel_dow[i][0].result(), flush=True)

      parallel_jobs = []
      results = []

      for index, row in csv_file.iterrows():
        pdb = row['Pdb']
        pdb_work = row['Pdb work']

        pdb = pdb.replace("+", "")
        pdb = pdb.replace("5000000000000.0", "5e12")

        pdb_work = pdb_work.replace("+", "")
        pdb_work = pdb_work.replace("5000000000000.0", "5e12")

        make_dir_pdb(work_dir, pdb)
        make_dir_pdb(work_dir, pdb_work)

        parallel_jobs.append(executor.submit(check_RMSD_result_algorithm,
                                             work_dir,
                                             row['Chains'],
                                             row['Match'],
                                             pdb,
                                             pdb_work,
                                             row['Changed chain'],
                                             row['Changed chain Pdb work']))

      for i in parallel_jobs:
        add_result = i.result()
        results.append(add_result)

      csv_file['RMSD'] = results
      csv_file.to_csv(file_name)


main()
