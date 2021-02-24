import os
import shutil
from ftplib import FTP

from general_utils.string_utils import get_float_between_ss, get_float_value
from general_utils.temp_utils import gen_dir, free_dir
from general_utils.terminal_utils import get_out


def get_all_emd_name():
  list_servers = [["ftp.wwpdb.org", "/pub/emdb/structures/"],
                  ["ftp.rcsb.org", "/pub/emdb/structures/"],
                  ["ftp.ebi.ac.uk", "/pub/databases/emdb/structures/"],
                  ["ftp.pdbj.org", "/pub/emdb/structures/"]]

  for i in list_servers:
    try:
      ftp = FTP()
      ftp.connect(i[0])
      ftp.login()
      ftp.cwd(i[1])
      files_list = ftp.nlst()

      result = []

      for filename in files_list:
        result.append(filename[4:])

      return result
    except Exception as e:
      pass
  return e


def get_mass_angstrom(map_path):
  map_real_path = os.path.abspath(map_path)
  path = gen_dir()
  # path = "./temp_map_mass"
  if os.path.exists(path):
    shutil.rmtree(path)
  os.mkdir(path)

  f = open(path + "/fit.cxc", "w+")
  f.write("volume #1 origin 0,0,0 \r\n")
  f.write("measure volume #1\r\n")
  f.write("exit")
  f.close()

  commands_real_path = os.path.abspath(path + "/fit.cxc")
  mass = 0
  error, exit_binary_text = get_out("chimerax", "--nogui", map_real_path, commands_real_path)
  if error != 0:
    free_dir(path)
    raise Exception("Error on try to get mass")
  text = exit_binary_text
  mass = get_float_between_ss(text, "Enclosed volume for surface (#1.1) =", "\n")
  free_dir(path)
  return mass


def get_mrc_level(map_path):
  map_real_path = os.path.abspath(map_path)
  path = gen_dir()
  # path = "./temp_map_mass"
  if os.path.exists(path):
    shutil.rmtree(path)
  os.mkdir(path)

  f = open(path + "/fit.cxc", "w+")
  f.write("exit")
  f.close()

  commands_real_path = os.path.abspath(path + "/fit.cxc")
  level = 0
  error, exit_binary_text = get_out("chimerax", "--nogui", map_real_path, commands_real_path)
  if error != 0:
    free_dir(path)
    raise Exception("Error on try to get mass")
  text = exit_binary_text
  level = get_float_between_ss(text, "at level", ",")
  free_dir(path)
  return level


def get_cube_len_angstrom(map_path):
  map_real_path = os.path.abspath(map_path)
  path = gen_dir()
  # path = "./temp_map_center"
  if os.path.exists(path):
    shutil.rmtree(path)

  os.mkdir(path)
  f = open(path + "/fit.cxc", "w+")
  f.write("volume #1 dumpHeader true \r\n")
  f.write("exit")
  f.close()

  commands_real_path = os.path.abspath(path + "/fit.cxc")

  _error, exit_binary_text = get_out("chimerax", "--nogui", map_real_path, commands_real_path)
  text = exit_binary_text
  x = get_float_value(text, 'xlen =', '\n')
  y = get_float_value(text, 'ylen =', '\n')
  z = get_float_value(text, 'zlen =', '\n')

  free_dir(path)

  return [x, y, z]
