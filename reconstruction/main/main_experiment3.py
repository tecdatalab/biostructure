import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")


def experiment_3():
  from experiment.experiment_3 import do_parallel_test_a

  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/data_experiment_3_a_v1".format(local_path), "result.csv", [3.5, 9.5])
  print("Finish")


if __name__ == '__main__':
  experiment_3()
