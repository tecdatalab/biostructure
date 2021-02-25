import os
import shutil
import tempfile

global_temp_dir = None


def free_dir(path_dir):
  shutil.rmtree(path_dir)


def gen_dir():
  if global_temp_dir == None:
    return os.path.abspath(tempfile.mkdtemp())

  if not os.path.exists(global_temp_dir):
    os.mkdir(global_temp_dir)

  return os.path.abspath(tempfile.mkdtemp(dir=global_temp_dir))


def clean_work_dir():
  if global_temp_dir != None:
    if os.path.exists(global_temp_dir):
      shutil.rmtree(global_temp_dir)
