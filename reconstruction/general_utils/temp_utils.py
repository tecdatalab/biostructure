import os
import shutil
import tempfile

global_temp_dir = None


def free_file(file_path):
  os.remove(file_path)

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


def gen_file_with_extension(extension):
  if global_temp_dir == None:
    _new_file, filename = tempfile.mkstemp()
    result = add_extension_in_file(filename, extension)
    return os.path.abspath(result)

  if not os.path.exists(global_temp_dir):
    try:
      os.mkdir(global_temp_dir)
    except:
      pass

  new_file, _filename = tempfile.mkstemp(dir=global_temp_dir)
  result = add_extension_in_file(_filename, extension)
  return os.path.abspath(result)


def clean_file(file):
  file_to_delete = open(file, 'w')
  file_to_delete.close()

def add_extension_in_file(original, extension):
  basename = os.path.basename(original)
  dirname = os.path.dirname(original)
  new_name = basename + extension
  new_path = os.path.join(dirname, new_name)
  os.rename(original, new_path)
  return new_path
