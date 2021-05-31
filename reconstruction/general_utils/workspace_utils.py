import os


def is_work_in_cluster():
  path = '/work/lcastillo/'
  isFile = os.path.isdir(path)
  return isFile
