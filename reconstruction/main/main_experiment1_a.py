import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import general_utils

from general_utils.workspace_utils import is_work_in_cluster
from csv_modules.csv_combine import combine_files_exp_1a
from general_utils.temp_utils import clean_work_dir

folder_work = "data_experiment_1_a_v1"


def experiment_1_a():
  from experiment.experiment_1_a import do_parallel_test_a
  if is_work_in_cluster():
    local_path = "/work/lcastillo"
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"

  print("Start")
  do_parallel_test_a("{0}/{1}".format(local_path, folder_work),
                     result_cvs_chain="result_chain.csv",
                     result_cvs_sequence="result_secuencial.csv",
                     result_cvs_struct="result_struct.csv",
                     resolution_range=[4, 6, 8, 10],
                     error_file="error_log_expe_1a.txt",
                     percentage_data_set=10,
                     file_checkpoint='check_expe_1a.pkl',
                     can_chain_test=3,
                     can_sequence_test=3,
                     can_struct_test=3,
                     add_to_ignore_files=False)
  print("Finish")


def union_test():
  if is_work_in_cluster():
    local_path = "/work/lcastillo/{0}".format(folder_work)
  else:
    local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/{0}".format(folder_work)

  combine_files_exp_1a('salida_struct.csv', 'salida_chain.csv', local_path)


if __name__ == '__main__':
  if is_work_in_cluster():
    general_utils.temp_utils.global_temp_dir = "/work/lcastillo/temp_exp_1a"
  else:
    general_utils.temp_utils.global_temp_dir = None
  clean_work_dir()
  experiment_1_a()
  # union_test()
