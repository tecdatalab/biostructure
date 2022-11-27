import os
import shutil

from general_utils.pdb_utils import move_pdb_center
from general_utils.string_utils import get_float_between_ss, get_float_value
from general_utils.temp_utils import gen_dir, free_dir
from general_utils.terminal_utils import get_out, execute_command
from general_utils.pdb_utils import align_pdb_file_1_in_2

def get_mass_angstrom(map_path, original):
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


def get_mrc_level(map_path,original,best_iterations_num=8):
  #best num is the number of times the sum and substraction is run
  #the bigger the number, the better
  #0.5 +-
  # get level from auxiliar function
  #add or substract 0.5 from that level to find which one is better
  #return the best number
  level=get_mrc_level_aux(map_path)
  if original!=None:
    level_tmp=level
    addition=0.5
    substraction=0.5
    #Generate simulated pdb with this level to calculate the metric con el pdb simulado(map_path) y con el original
    #check if the metric is correct if its ok, then return the result
    low_actual_value=level_tmp
    high_actual_value=level_tmp
    level_to_return=level
    #get the current metric for the pdbs
    pdb_simulated = get_mrc_to_pdb_aux(level, map_path, "simulated_1yfq.pdb")
    last_metric=align_pdb_file_1_in_2(pdb_simulated,original).RMSDAfterRefinement #initiates the last metric variable
    if os.path.isfile("/tmp/simulated_1yfq.pdb"):
      os.remove("/tmp/simulated_1yfq.pdb")
    for i in range(best_iterations_num):
      if os.path.isfile("/tmp/low_1yfq.pdb"):
        os.remove("/tmp/low_1yfq.pdb")
      if os.path.isfile("/tmp/high_1yfq.pdb"):
        os.remove("/tmp/high_1yfq.pdb")
      if metric_low!=float("inf"):
        pdb_simulated_low=get_mrc_to_pdb_aux(low_actual_value,map_path,"low_1yfq.pdb")
      if metric_high!=float("inf"):
        pdb_simulated_high=get_mrc_to_pdb_aux(high_actual_value,map_path,"high_1yfq.pdb")
      try:
        if metric_low==float("inf"):
          raise Exception
        metric_low=align_pdb_file_1_in_2(pdb_simulated_low,original).RMSDAfterRefinement
      except:
        metric_low=float("inf")
        substraction=0
      try:
        if metric_high==float("inf"):
          raise Exception
        metric_high=align_pdb_file_1_in_2(pdb_simulated_high,original).RMSDAfterRefinement
      except:
        metric_high=float("inf")
        addition=0
      if metric_low <= 0 or metric_high<=0:
        return min([metric_low,metric_high])
      elif metric_low==float("inf") and metric_high==float("inf"):
        break
      elif metric_low<=metric_high:
        if metric_low<=last_metric:
          last_metric=metric_low
          level_to_return=low_actual_value
      else:
        if metric_high<=last_metric:
          last_metric=metric_high
          level_to_return=high_actual_value

      low_actual_value-=substraction
      high_actual_value+=addition

    level=level_to_return
    if os.path.isfile("/tmp/low_1yfq.pdb"):
      os.remove("/tmp/low_1yfq.pdb")
    if os.path.isfile("/tmp/high_1yfq.pdb"):
      os.remove("/tmp/high_1yfq.pdb")
  return level

def get_mrc_to_pdb_aux(level, mrc_path,name_file):
  pdb_path="/tmp/"+name_file
  execute_command(
    "../binaries/MAINMAST/MainmastC -i {0} -t {1} -r 50 > {2}".format(
      mrc_path, level, pdb_path))

  from general_utils.pdb_utils import only_first_model
  only_first_model(pdb_path)


  move_pdb_center(pdb_path)

  return pdb_path


def get_mrc_level_aux(map_path):
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


def mrc_to_pdb_test(mrc_path,pdb_result_path,original_path=None, threshold_mrc=0.0, clean=False):
  level = get_mrc_level(mrc_path, original_path)

  if (threshold_mrc == 0 or threshold_mrc == 0.0):
    threshold_mrc = level

  execute_command(
    "../binaries/MAINMAST/MainmastC -i {0} -t {1} -r 50 > {2}".format(
      mrc_path, threshold_mrc, pdb_result_path))
  # free_dir(temp_dir)

  if clean:
    from general_utils.pdb_utils import only_first_model
    only_first_model(pdb_result_path)
  move_pdb_center(pdb_result_path)


def mrc_to_pdb(mrc_path, pdb_result_path, threshold_mrc=0.0, clean=False):
  # temp_dir = gen_dir()

  #execute_command("echo 1|../binaries/Situs_3.1/bin/map2map {0} {1}/tempPDB.situs".format(mrc_path, temp_dir))

<<<<<<< HEAD
  level = get_mrc_level(mrc_path, original_pdb)
=======
  #Add flag, if original archive is here, find the best level, else, dont find best level
  level = get_mrc_level(mrc_path)
>>>>>>> em-reconstruction

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
