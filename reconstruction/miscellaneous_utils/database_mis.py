import os

import gridfs
import pickle

# Constans
from general_utils.temp_utils import gen_dir, free_dir, gen_file

ZIPJSON_KEY = 'base64(zip(o))'


# Big Files
def put_bigfile_db(db, data_save):
  file_write = gen_file()
  output = open(file_write, 'wb')
  pickle.dump(data_save, output)
  output.close()

  output = open(file_write, 'rb')
  data_stream = output.read()

  fs = gridfs.GridFS(db)
  fid = fs.put(data_stream)

  output.close()
  os.remove(file_write)

  return fid


def get_bigfile_db(db, key_file):
  fs = gridfs.GridFS(db)
  fs.exists(key_file)
  data_stream = fs.get(key_file).read()

  file_write = gen_file()
  output = open(file_write, 'wb')
  pickle.dump(data_stream, output)
  output.close()

  output = open(file_write, 'rb')
  outputdata = pickle.load(output)

  output.close()
  os.remove(file_write)

  return outputdata


# Graph


# Files
def get_dicc_chains_files_info(pdb_id, path=None):
  from general_utils.pdb_utils import make_pdb_dir_temp

  if path == None:
    path = gen_dir()

  make_pdb_dir_temp(path, pdb_id)
  dicc_add = {}
  path_pdb = os.path.join(path, pdb_id)
  for chain_file in os.listdir(path_pdb):
    path_chain = os.path.join(path_pdb, chain_file)
    with open(path_chain) as f:
      name_chain = path_chain.split(".")[0]
      name_chain = name_chain.split("_")[1]
      dicc_add[name_chain] = "".join(f.readlines())

  free_dir(path)
  return dicc_add
