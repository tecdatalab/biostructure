import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import general_utils

from general_utils.workspace_utils import is_work_in_cluster
from general_utils.temp_utils import clean_work_dir
from csv_modules.csv_combine import combine_files_exp_1

folder_work = "data_experiment_1_d_v1"


def experiment_1_d():
  from experiment.experiment_1_d import do_parallel_test_a
  if is_work_in_cluster():
    local_path = "/work/lcastillo"
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"

  print("Start")
  do_parallel_test_a("{0}/{1}".format(local_path, folder_work),
                     "result.csv",
                     [4, 6, 8, 10],
                     #pdbs_work=['2ian', '4fmi', '6l7o', '2nx5', '6ytk', '2df7', '6gej', '4ind', '3u8k', '6m6h', '2qjh', '6rdm', '6z86', '5gip', '3glc', '6xky', '3r8r', '6w09', '5no4', '5w66'],
                     pdbs_work=['2ian'],
                     error_file="error_log_expe_1_d.txt")
  print("Finish")


def union_test():
  if is_work_in_cluster():
    local_path = "/work/lcastillo/{0}".format(folder_work)
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/{0}".format(folder_work)

  combine_files_exp_1('salida_exper1_d_v1.csv', local_path)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_exp_1"
  else:
    general_utils.temp_utils.global_temp_dir = None
  clean_work_dir()
  experiment_1_d()
  # result = get_similar_pdb("6LFS")
  # union_test()
