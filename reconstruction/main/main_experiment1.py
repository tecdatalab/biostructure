import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import general_utils
from general_utils.temp_utils import clean_work_dir
from csv_modules.csv_combine import combine_files_exp_1

folder_work = "data_experiment_1_v2"

def experiment_1():
  from experiment.experiment_1 import do_parallel_test_a
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/{1}".format(local_path, folder_work), "result.csv", [4, 6, 8, 10],
                     range_incompleteness=[10.0, 15.0], can_try_experiments=10, add_to_ignore_files=False,
                     error_file="error_log_expe_1.txt")
  print("Finish")


def union_test():
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/{0}".format(folder_work)
  # local_path = "/work/lcastillo/{0}".format(folder_work)
  combine_files_exp_1('salida_exper1_v2.csv',
                      local_path)


if __name__ == '__main__':
  general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_exp_1"
  general_utils.temp_utils.global_temp_dir = None
  clean_work_dir()
  #experiment_1()
  #result = get_similar_pdb("6LFS")
  union_test()
