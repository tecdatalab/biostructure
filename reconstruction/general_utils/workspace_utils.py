import os


def is_work_in_cluster():
  from constans.directory_constans import CLUSTER_WORK_PATH
  path = CLUSTER_WORK_PATH
  isFile = os.path.isdir(path)
  return isFile
