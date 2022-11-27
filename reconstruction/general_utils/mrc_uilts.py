import os
import shutil

from general_utils.pdb_utils import move_pdb_center
from general_utils.string_utils import get_float_between_ss, get_float_value
from general_utils.temp_utils import gen_dir, free_dir
from general_utils.terminal_utils import get_out, execute_command


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


def mrc_to_pdb(mrc_path, pdb_result_path, threshold_mrc=0.0, clean=False):
  # temp_dir = gen_dir()

  #execute_command("echo 1|../binaries/Situs_3.1/bin/map2map {0} {1}/tempPDB.situs".format(mrc_path, temp_dir))

  level = get_mrc_level(mrc_path, original_pdb)

  if(threshold_mrc==0 or threshold_mrc==0.0):
    threshold_mrc = level

  execute_command(
    "../binaries/MAINMAST/MainmastC -i {0} -t {1} -r 50 > {2}".format(
      mrc_path, threshold_mrc, pdb_result_path))
  # free_dir(temp_dir)

  if clean:
    from general_utils.pdb_utils import only_first_model
    only_first_model(pdb_result_path)
  move_pdb_center(pdb_result_path)
