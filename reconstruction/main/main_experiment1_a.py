


def experiment_1_a():
  from experiment.experiment_1_a import do_parallel_test_a

  local_path = "/home/lcastillo98/Documents/git_projects/biostructure/reconstruction"
  # local_path = "/work/lcastillo"
  print("Start")
  do_parallel_test_a("{0}/data_experiment_1_a_a_v1".format(local_path), result_cvs_chain="result_chain.csv",
                     result_cvs_struct="result_struct.csv", resolution_range=[3.5, 9.5],
                     error_file="error_log_expe_1a.txt", percentage_data_set=10, file_checkpoint='check_expe_1a.pkl',
                     can_chain_test=3, can_struct_test=3)
  print("Finish")


if __name__ == '__main__':
  experiment_1_a()
