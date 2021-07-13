import sys
import pathlib

from experiment.utils_general import check_RMSD_result_algorithm, make_dir_pdb

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import pandas as pd
from ast import literal_eval
import os

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

      parallel_jobs = []
      results = []

      for index, row in csv_file.iterrows():

        make_dir_pdb(work_dir, row['Pdb'])
        make_dir_pdb(work_dir, row['Pdb work'])

        parallel_jobs.append(executor.submit(check_RMSD_result_algorithm,
                                             work_dir,
                                             row['Chains'],
                                             row['Match'],
                                             row['Pdb'],
                                             row['Pdb work'],
                                             row['Changed chain'],
                                             row['Changed chain Pdb work']))

      for i in parallel_jobs:
        add_result = i.result()
        results.append(add_result)

      csv_file['RMSD'] = results
      csv_file.to_csv(file_name)

main()
