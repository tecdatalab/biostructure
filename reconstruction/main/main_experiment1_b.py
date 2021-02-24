import sys
import pathlib


sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from general_utils.temp_utils import clean_work_dir
from csv_modules.csv_combine import combine_files_exp_1a
folder_work = "data_experiment_1_b_v1"


def experiment_1_b():
  from experiment.experiment_1_b import do_parallel_test
  clean_work_dir()
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test("{0}/{1}".format(local_path, folder_work),
                     result_cvs_chain="result_chain.csv",
                     result_cvs_struct="result_struct.csv",
                     result_cvs_secuencial="result_secuencial.csv",
                     resolution_range=[3.5, 9.5],
                     error_file="error_log_expe_1b.txt", percentage_data_set=10, file_checkpoint='check_expe_1b.pkl',
                     can_chain_test=3, can_struct_test=3, can_secuencial_test=3, add_to_ignore_files=False)
  print("Finish")


def union_test():
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/{0}".format(folder_work)
  # local_path = "/work/lcastillo/{0}".format(folder_work)
  combine_files_exp_1a('salida_struct.csv', 'salida_chain.csv', local_path)


if __name__ == '__main__':
  experiment_1_b()
  # union_test()
