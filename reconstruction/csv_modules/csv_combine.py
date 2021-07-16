import os
import pandas as pd
from os import listdir
from ast import literal_eval
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from experiment.utils_general import make_dir_pdb, check_RMSD_result_algorithm
from general_utils.workspace_utils import is_work_in_cluster


def find_csv_filenames(path_to_dir, result, suffix=".csv"):
  filenames = listdir(path_to_dir)
  for filename in filenames:
    if os.path.isdir('{}/{}'.format(path_to_dir, filename)):
      find_csv_filenames('{}/{}'.format(path_to_dir, filename), result, suffix)
    else:
      if filename.endswith(suffix):
        result.append('{}/{}'.format(path_to_dir, filename))


def count_no_match_chains(x):
  can_chains_f = x['number_chains']
  can_chains_t = x['number_test_chains']
  can_match = x['total_matched']

  return min(can_chains_f, can_chains_t) - can_match


def count_error(x):
  cont = 0
  for i in x:
    if i[0] != i[1]:
      cont += 1
  return cont


def count_ok(x):
  cont = 0
  for i in x:
    if i[0] == i[1]:
      cont += 1
  return cont

def percentage_wrong(x):
  can_chains = x['total_matched']
  can_wrong = x['number_wrong_chains']

  if can_chains == 0:
    return 1.0
  return can_wrong/can_chains

def percentage_ok(x):
  can_chains = x['total_matched']
  can_ok = x['number_ok_chains']

  if can_chains == 0:
    return 1.0
  return can_ok/can_chains

def percentage_no_match(x):
  can_chains_f = x['number_chains']
  can_chains_t = x['number_test_chains']
  can_no_match_chains= x['number_no_match_chains']

  if min(can_chains_f, can_chains_t) == 0:
    return 1.0
  return can_no_match_chains/min(can_chains_f, can_chains_t)

def percentage_match(x):
  can_match = x['total_matched']
  can_chains_f = x['number_chains']
  can_chains_t = x['number_test_chains']

  if min(can_chains_f, can_chains_t) == 0:
    return 1.0
  return can_match/min(can_chains_f, can_chains_t)


def combine_files_exp_1(exit_file, parent_path):
  result = []
  find_csv_filenames(parent_path, result)
  # print(result)

  # combine all files in the list
  combined_csv = pd.concat([pd.read_csv(f, converters={"Father Chains": literal_eval,
                                                       "Test Chains": literal_eval,
                                                       "Match": literal_eval}) for f in result])

  # Process data
  combined_csv['number_chains'] = combined_csv['Father Chains'].apply(lambda x: len(x))
  combined_csv['number_test_chains'] = combined_csv['Test Chains'].apply(lambda x: len(x))
  combined_csv['total_matched'] = combined_csv['Match'].apply(lambda x: len(x))
  combined_csv['number_wrong_chains'] = combined_csv['Match'].apply(count_error)
  combined_csv['number_ok_chains'] = combined_csv['Match'].apply(count_ok)
  combined_csv['number_no_match_chains'] = combined_csv.apply(count_no_match_chains, axis=1)
  combined_csv['percentage_wrong'] = combined_csv.apply(percentage_wrong, axis=1)
  combined_csv['percentage_ok'] = combined_csv.apply(percentage_ok, axis=1)
  combined_csv['percentage_no_match'] = combined_csv.apply(percentage_no_match, axis=1)
  combined_csv['percentage_match'] = combined_csv.apply(percentage_match, axis=1)

  # export to csv
  combined_csv.to_csv(exit_file, index=False, encoding='utf-8-sig')


def combine_files_exp_1a(exit_file_struct, exit_file_chain, parent_path):
  result_struct = []
  result_chain = []
  find_csv_filenames(parent_path, result_struct, suffix="struct.csv")
  find_csv_filenames(parent_path, result_chain, suffix="chain.csv")
  # print(result_struct, result_chain)

  # combine all files in the list
  combined_csv_struct = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval,
                                                       "Work Chains": literal_eval}) for f in result_struct])
  combined_csv_struct['number_chains'] = combined_csv_struct['Chains'].apply(lambda x: len(x))
  combined_csv_struct['number_work_chains'] = combined_csv_struct['Work Chains'].apply(lambda x: len(x))

  combined_csv_chains = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval}) for f in result_chain])
  combined_csv_chains['number_chains'] = combined_csv_chains['Chains'].apply(lambda x: len(x))

  # export to csv
  combined_csv_struct.to_csv(exit_file_struct, index=False, encoding='utf-8-sig')
  combined_csv_chains.to_csv(exit_file_chain, index=False, encoding='utf-8-sig')


def add_RSM(dataFrame, executor, work_dir):
  col_names_dowloand = dataFrame['Pdb'].tolist()
  col_names_dowloand += dataFrame['Pdb work'].tolist()
  col_names_dowloand = np.unique(col_names_dowloand).tolist()

  parallel_dow = []
  for i in col_names_dowloand:
    parallel_dow.append([executor.submit(make_dir_pdb,
                                         work_dir,
                                         i),
                         i])

    for i in range(len(parallel_dow)):
      print(parallel_dow[i][1], i / len(parallel_dow), flush=True)
      print(parallel_dow[i][0].result(), flush=True)

    parallel_jobs = []
    results = []

    for index, row in dataFrame.iterrows():
      pdb = row['Pdb']
      pdb_work = row['Pdb work']

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

    dataFrame['RMSD'] = results
    return dataFrame

def combine_files_exp_1b(exit_file_struct, exit_file_chain, exit_file_secuencial, parent_path):
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      result_struct = []
      result_chain = []
      result_sequence = []
      find_csv_filenames(parent_path, result_struct, suffix="struct.csv")
      find_csv_filenames(parent_path, result_chain, suffix="chain.csv")
      find_csv_filenames(parent_path, result_sequence, suffix="secuencial.csv")
      # print(result_struct, result_chain)

      # combine all files in the list
      combined_csv_struct = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval,
                                                                  "Work Chains": literal_eval,
                                                                  "Match": literal_eval}) for f in result_struct])

      combined_csv_sequence = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval,
                                                                  "Match": literal_eval}) for f in result_sequence])

      combined_csv_chain = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval,
                                                                  "Match": literal_eval}) for f in result_chain])



      combined_csv_struct['number_chains'] = combined_csv_struct['Chains'].apply(lambda x: len(x))
      combined_csv_struct['number_test_chains'] = combined_csv_struct['Work Chains'].apply(lambda x: len(x))

      combined_csv_sequence['number_chains'] = combined_csv_sequence['Chains'].apply(lambda x: len(x))
      combined_csv_sequence['number_test_chains'] = combined_csv_sequence['Chains'].apply(lambda x: len(x))
      combined_csv_sequence['total_matched'] = combined_csv_sequence['Match'].apply(lambda x: len(x))
      combined_csv_sequence['number_wrong_chains'] = combined_csv_sequence['Match'].apply(count_error)
      combined_csv_sequence['number_ok_chains'] = combined_csv_sequence['Match'].apply(count_ok)
      combined_csv_sequence['number_no_match_chains'] = combined_csv_sequence.apply(count_no_match_chains, axis=1)
      combined_csv_sequence['percentage_wrong'] = combined_csv_sequence.apply(percentage_wrong, axis=1)
      combined_csv_sequence['percentage_ok'] = combined_csv_sequence.apply(percentage_ok, axis=1)
      combined_csv_sequence['percentage_no_match'] = combined_csv_sequence.apply(percentage_no_match, axis=1)
      combined_csv_sequence['percentage_match'] = combined_csv_sequence.apply(percentage_match, axis=1)

      combined_csv_chain['number_chains'] = combined_csv_chain['Chains'].apply(lambda x: len(x))
      combined_csv_chain['number_test_chains'] = combined_csv_chain['Chains'].apply(lambda x: len(x))
      combined_csv_chain['total_matched'] = combined_csv_chain['Match'].apply(lambda x: len(x))
      combined_csv_chain['number_wrong_chains'] = combined_csv_chain['Match'].apply(count_error)
      combined_csv_chain['number_ok_chains'] = combined_csv_chain['Match'].apply(count_ok)
      combined_csv_chain['number_no_match_chains'] = combined_csv_chain.apply(count_no_match_chains, axis=1)
      combined_csv_chain['percentage_wrong'] = combined_csv_chain.apply(percentage_wrong, axis=1)
      combined_csv_chain['percentage_ok'] = combined_csv_chain.apply(percentage_ok, axis=1)
      combined_csv_chain['percentage_no_match'] = combined_csv_chain.apply(percentage_no_match, axis=1)
      combined_csv_chain['percentage_match'] = combined_csv_chain.apply(percentage_match, axis=1)

      combined_csv_struct['number_chains'] = combined_csv_struct['Chains'].apply(lambda x: len(x))
      combined_csv_struct['number_test_chains'] = combined_csv_struct['Work Chains'].apply(lambda x: len(x))
      combined_csv_struct['total_matched'] = combined_csv_struct['Match'].apply(lambda x: len(x))
      combined_csv_struct['percentage_match'] = combined_csv_struct.apply(percentage_match, axis=1)

      if is_work_in_cluster():
        add_RSM(combined_csv_sequence, executor, "/work/lcastillo/RMSD")
      else:
        add_RSM(combined_csv_sequence, executor, "./RMSD")

      # export to csv
      combined_csv_struct.to_csv(exit_file_struct, index=False, encoding='utf-8-sig')
      combined_csv_chain.to_csv(exit_file_chain, index=False, encoding='utf-8-sig')
      combined_csv_sequence.to_csv(exit_file_secuencial, index=False, encoding='utf-8-sig')



def combine_files_exp_1c(exit_file, parent_path):
  result = []
  find_csv_filenames(parent_path, result)
  # print(result)

  # combine all files in the list
  combined_csv = pd.concat([pd.read_csv(f, converters={"Chains": literal_eval,
                                                       "Nodes father": literal_eval,
                                                       "Nodes test": literal_eval,
                                                       "Match": literal_eval}) for f in result])

  # Process data
  combined_csv['number_chains'] = combined_csv['Nodes father'].apply(lambda x: len(x))
  combined_csv['number_test_chains'] = combined_csv['Nodes test'].apply(lambda x: len(x))
  combined_csv['total_matched'] = combined_csv['Match'].apply(lambda x: len(x))
  combined_csv['number_wrong_chains'] = combined_csv['Match'].apply(count_error)
  combined_csv['number_ok_chains'] = combined_csv['Match'].apply(count_ok)
  combined_csv['number_no_match_chains'] = combined_csv.apply(count_no_match_chains, axis=1)
  combined_csv['percentage_wrong'] = combined_csv.apply(percentage_wrong, axis=1)
  combined_csv['percentage_ok'] = combined_csv.apply(percentage_ok, axis=1)
  combined_csv['percentage_no_match'] = combined_csv.apply(percentage_no_match, axis=1)
  combined_csv['percentage_match'] = combined_csv.apply(percentage_match, axis=1)


  # export to csv
  combined_csv.to_csv(exit_file, index=False, encoding='utf-8-sig')
