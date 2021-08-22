import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
from general_utils.database_utils import get_chains_pdb_db
import time
from experiment.utils_general import remove_get_dirs, check_RMSD_result_algorithm, check_RMSD_result_all
import general_utils
from general_utils.workspace_utils import is_work_in_cluster
from general_utils.temp_utils import clean_work_dir, gen_dir, free_dir
import pandas as pd
from ast import literal_eval
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import os

# file_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/data_experiment_1_b_v1_exe_1/3sdd/result_struct.csv"
file_path = "/home/lcastillo/workspaces/project_biostructure/salida_struct_exe_3.csv"
# temp_csv = "./temp_struct.csv"

temp_csv = "/home/lcastillo/workspaces/project_biostructure/salida_struct_exe_3_RMSD_repare.csv"
temp_csv_v1 = "/home/lcastillo/workspaces/project_biostructure/salida_struct_exe_3_RMSD_repare_v1.csv"

RMSD_dir = "/work/lcastillo/chain_pdb/"

def experiment_1_b_RMSD_chain():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      if not os.path.exists(temp_csv):
        csv_file = pd.read_csv(file_path, converters={"Chains": literal_eval,
                                                      "Match": literal_eval,
                                                      "Pdb": str,
                                                      "Pdb work": str,
                                                      "Changed chain Pdb work": str,
                                                      "Changed chain": str,
                                                      "RMSD": literal_eval
                                                      })

        for index, row in csv_file.iterrows():
          csv_file.loc[index, 'RMSD'] = False

        csv_file.to_csv(temp_csv, index=False)

      csv_file = pd.read_csv(temp_csv, converters={"Chains": literal_eval,
                                                   "Match": literal_eval,
                                                   "Pdb": str,
                                                   "Pdb work": str,
                                                   "Changed chain Pdb work": str,
                                                   "Changed chain": str,
                                                   "RMSD": literal_eval
                                                   })

      new_RMSD_list_update = []
      list_update = []
      con_update = 0
      for row in csv_file.iterrows():

        if row[1].RMSD == False:
          list_update.append([True, executor.submit(calculate_RMSD_chain, row)])
          new_RMSD_list_update.append(False)
          con_update += 1
        else:
          list_update.append([False])
          new_RMSD_list_update.append(row[1].RMSD)

      total_updated = 0

      while total_updated < con_update:
        time.sleep(30)
        pos = 0
        for _row in csv_file.iterrows():
          if list_update[pos][0] and list_update[pos][1].done():
            try:
              RMSD_new = list_update[pos][1].result()
            except:
              RMSD_new = False
            new_RMSD_list_update[pos] = RMSD_new
            total_updated += 1
            print("To repare:", total_updated / con_update, total_updated, con_update, flush=True)

            list_update[pos][0] = False

          pos += 1

        if total_updated != 0 and total_updated % 10 == 0:
          csv_file["RMSD"] = new_RMSD_list_update
          csv_file.to_csv(temp_csv_v1, index=False)
          os.remove(temp_csv)
          os.rename(temp_csv_v1, temp_csv)

      csv_file["RMSD"] = new_RMSD_list_update
      csv_file.to_csv(temp_csv_v1, index=False)
      os.remove(temp_csv)
      os.rename(temp_csv_v1, temp_csv)

      print("Finish")


def experiment_1_b_RMSD_struct():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      if not os.path.exists(temp_csv):
        csv_file = pd.read_csv(file_path, converters={"Chains": literal_eval,
                                                      "Work Chains": literal_eval,
                                                      "Match": literal_eval,
                                                      "Pdb": str,
                                                      "Pdb work": str,
                                                      "RMSD": literal_eval
                                                      })

        for index, row in csv_file.iterrows():
          csv_file.loc[index, 'RMSD'] = False

        csv_file.to_csv(temp_csv, index=False)

      csv_file = pd.read_csv(temp_csv, converters={"Chains": literal_eval,
                                                   "Work Chains": literal_eval,
                                                   "Match": literal_eval,
                                                   "Pdb": str,
                                                   "Pdb work": str,
                                                   "RMSD": literal_eval
                                                   })

      new_RMSD_list_update = []
      list_update = []
      con_update = 0
      for row in csv_file.iterrows():

        if row[1].RMSD == False:
          list_update.append([True, executor.submit(calculate_RMSD_struct, row)])
          new_RMSD_list_update.append(False)
          con_update += 1
        else:
          list_update.append([False])
          new_RMSD_list_update.append(row[1].RMSD)

      total_updated = 0

      while total_updated < con_update:
        time.sleep(30)
        pos = 0
        for _row in csv_file.iterrows():
          if list_update[pos][0] and list_update[pos][1].done():
            try:
              RMSD_new = list_update[pos][1].result()
            except:
              RMSD_new = False
            new_RMSD_list_update[pos] = RMSD_new
            total_updated += 1
            print("To repare:", total_updated / con_update, total_updated, con_update, flush=True)

            list_update[pos][0] = False

          pos += 1

        if total_updated != 0 and total_updated % 10 == 0:
          csv_file["RMSD"] = new_RMSD_list_update
          csv_file.to_csv(temp_csv_v1, index=False)
          os.remove(temp_csv)
          os.rename(temp_csv_v1, temp_csv)

      csv_file["RMSD"] = new_RMSD_list_update
      csv_file.to_csv(temp_csv_v1, index=False)
      os.remove(temp_csv)
      os.rename(temp_csv_v1, temp_csv)

      print("Finish")


def calculate_RMSD_chain(row):
  result = check_RMSD_result_algorithm(RMSD_dir,
                                       get_chains_pdb_db(row[1][0]),
                                       row[1].Match,
                                       row[1][0],
                                       row[1][1],
                                       row[1][14],
                                       row[1][15])

  return result


def calculate_RMSD_struct(row):
  result = check_RMSD_result_all(RMSD_dir,
                                 get_chains_pdb_db(row[1][0]),
                                 get_chains_pdb_db(row[1][1]),
                                 row[1][11],
                                 row[1][0],
                                 row[1][1])

  return result


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_exp_1b_RMSD"
  else:
    general_utils.temp_utils.global_temp_dir = None

  clean_work_dir()
  # experiment_1_b_RMSD_chain()
  experiment_1_b_RMSD_struct()
