import os
import pandas as pd
from os import listdir
from ast import literal_eval
import numpy as np

def find_csv_filenames(path_to_dir, result, suffix=".csv"):
  filenames = listdir(path_to_dir)
  for filename in filenames:
    if os.path.isdir('{}/{}'.format(path_to_dir, filename)):
      find_csv_filenames('{}/{}'.format(path_to_dir, filename), result, suffix)
    else:
      if filename.endswith(suffix):
        result.append('{}/{}'.format(path_to_dir, filename))


def count_no_match_chains(x):
  f_chains = x['Father Chains']
  t_chains = x['Test Chains']
  match = x['Match']

  dicc_chains = {}
  test = []
  check_test = []
  for i in range(len(f_chains)):
    dicc_chains[f_chains[i]] = i+1

  for i in t_chains:
    test.append(dicc_chains[i])

  for i in match:
    check_test.append(i[1])

  result = np.setdiff1d(test, check_test)

  return result.shape[0]


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
  can_chains = x['number_test_chains']
  can_wrong = x['number_wrong_chains']

  return can_wrong/can_chains

def percentage_ok(x):
  can_chains = x['number_test_chains']
  can_ok = x['number_ok_chains']

  return can_ok/can_chains

def percentage_no_match(x):
  can_chains = x['number_test_chains']
  can_no_match_chains= x['number_no_match_chains']

  return can_no_match_chains/can_chains


def combine_files(exit_file, parent_path):
  result = []
  find_csv_filenames(parent_path, result)
  print(result)

  # combine all files in the list
  combined_csv = pd.concat([pd.read_csv(f, converters={"Father Chains": literal_eval,
                                                       "Test Chains": literal_eval,
                                                       "Match": literal_eval}) for f in result])

  # Process data
  combined_csv['number_chains'] = combined_csv['Father Chains'].apply(lambda x: len(x))
  combined_csv['number_test_chains'] = combined_csv['Test Chains'].apply(lambda x: len(x))
  combined_csv['number_wrong_chains'] = combined_csv['Match'].apply(count_error)
  combined_csv['number_ok_chains'] = combined_csv['Match'].apply(count_ok)
  combined_csv['number_no_match_chains'] = combined_csv.apply(count_no_match_chains, axis=1)
  combined_csv['percentage_wrong'] = combined_csv.apply(percentage_wrong, axis=1)
  combined_csv['percentage_ok'] = combined_csv.apply(percentage_ok, axis=1)
  combined_csv['percentage_no_match'] = combined_csv.apply(percentage_no_match, axis=1)

  # export to csv
  combined_csv.to_csv(exit_file, index=False, encoding='utf-8-sig')
