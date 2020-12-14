import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from experiment.utils_experiment_1 import get_similar_pdb
from csv_modules.csv_combine import combine_files


def experiment_1():
  from experiment.experiment_1 import do_parallel_test_a, do_parallel_test_b

  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/data_experiment_1_a_v2".format(local_path), "result.csv", [3.5, 9.5],
                     range_incompleteness=[10.0, 15.0], can_try_experiments=10)
  print("Finish")


def union_test():
  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction/data_experiment_1_a"
  # local_path = "/work/lcastillo/data_experiment_1_a"
  combine_files('salida.csv',
                local_path)


if __name__ == '__main__':
  experiment_1()
  #result = get_similar_pdb("6LFS")
