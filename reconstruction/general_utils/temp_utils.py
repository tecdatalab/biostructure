import os
import shutil
import tempfile

global_temp_dir = None


def free_dir(path_dir):
  try:
    shutil.rmtree(path_dir)
  except:
    pass


def gen_dir():
  if global_temp_dir == None:
    return os.path.abspath(tempfile.mkdtemp())

  if not os.path.exists(global_temp_dir):
    try:
      os.mkdir(global_temp_dir)
    except:
      pass


  return os.path.abspath(tempfile.mkdtemp(dir=global_temp_dir))


def clean_work_dir():
  if global_temp_dir != None:
    if os.path.exists(global_temp_dir):
      try:
        shutil.rmtree(global_temp_dir)
      except:
        pass

      try:
        os.mkdir(global_temp_dir)
      except:
        pass


def gen_file():
  if global_temp_dir == None:
    _new_file, filename = tempfile.mkstemp()
    return os.path.abspath(filename)

  if not os.path.exists(global_temp_dir):
    try:
      os.mkdir(global_temp_dir)
    except:
      pass

  new_file, _filename = tempfile.mkstemp(dir=global_temp_dir)
  return os.path.abspath(_filename)
