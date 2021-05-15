import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import general_utils
from general_utils.workspace_utils import is_work_in_cluster
from general_utils.temp_utils import clean_work_dir
from csv_modules.csv_combine import combine_files_exp_1a

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--number_exe', help='foo help')
args = parser.parse_args()

if args.number_exe is None:
  raise ("Can not run")


folder_work = "data_experiment_1_c_v1_exe_{}".format(args.number_exe)


def experiment_1_c():
  from experiment.experiment_1_c import do_parallel_test
  if is_work_in_cluster():
    local_path = "/work/lcastillo"
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"

  print("Start")
  do_parallel_test("{0}/{1}".format(local_path, folder_work),
                   result_cvs="result.csv",
                   resolution_range=[4, 6, 8, 10],
                   error_file="error_log_expe_1c_exe_{}.txt".format(args.number_exe), percentage_data_set=100,
                   file_checkpoint='check_expe_1c.pkl',
                   add_to_ignore_files=False, can_groups=3, min_can_chains=6)
  print("Finish")


def union_test():
  if is_work_in_cluster():
    local_path = "/work/lcastillo/{0}".format(folder_work)
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/{0}".format(folder_work)

  combine_files_exp_1a('salida_struct.csv', 'salida_chain.csv', local_path)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_exp_1c_exe_{}".format(args.number_exe)
  else:
    general_utils.temp_utils.global_temp_dir = None

  clean_work_dir()
  experiment_1_c()
  # union_test()
