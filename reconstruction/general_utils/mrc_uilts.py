import os
import shutil

from general_utils.pdb_utils import move_pdb_center
from general_utils.string_utils import get_float_between_ss, get_float_value
from general_utils.temp_utils import gen_dir, free_dir, gen_file, free_file, gen_file_with_extension, clean_file
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


def get_mrc_level(map_path, original, best_iterations_num=8, break_sim=0.35, addition_value=0.5,
                  substraction_value=0.5):
  # best num is the number of times the sum and substraction is run
  # the bigger the number, the better
  # 0.5 +-
  # get level from auxiliar function
  # add or substract 0.5 from that level to find which one is better
  # return the best number
  level = get_mrc_level_aux(map_path)
  if original != None:
    addition = addition_value
    substraction = substraction_value
    # Generate simulated pdb with this level to calculate the metric con el pdb simulado(map_path) y con el original
    # check if the metric is correct if its ok, then return the result
    low_actual_value = level - substraction
    high_actual_value = level + addition
    level_to_return = level
    # get the current metric for the pdbs

    # Create file for simulate pdb
    simulate_path = gen_file_with_extension(".pdb")
    low_simulate_path = gen_file_with_extension(".pdb")
    high_simulate_path = gen_file_with_extension(".pdb")
    print(simulate_path)
    print(low_simulate_path)
    print(high_simulate_path)

    get_mrc_to_pdb_aux(level, map_path, simulate_path)
    aling_result = align_pdb_file_1_in_2(simulate_path, original)
    print("Print 1", level)
    print(aling_result.RMSDAfterRefinement)
    print(aling_result.PercentageAtomsAlignedAfterRefinement)
    print(aling_result.RawAlignmentScore)
    print(aling_result.NumberAlignedAtomsAfterRefinement)
    last_metric = aling_result.RMSDAfterRefinement  # initiates the last metric variable
    metric_low = last_metric
    metric_high = last_metric

    if last_metric <= break_sim:
      return level

    free_file(simulate_path)
    i = 0
    while i < best_iterations_num:
      if metric_low != float("inf"):
        pdb_simulated_low = get_mrc_to_pdb_aux(low_actual_value, map_path, low_simulate_path)
      if metric_high != float("inf"):
        pdb_simulated_high = get_mrc_to_pdb_aux(high_actual_value, map_path, high_simulate_path)

      try:
        if metric_low == float("inf"):
          raise Exception
        aling_result = align_pdb_file_1_in_2(pdb_simulated_low, original)

        print("Print 2", low_actual_value)
        print(aling_result.RMSDAfterRefinement)
        print(aling_result.PercentageAtomsAlignedAfterRefinement)
        print(aling_result.RawAlignmentScore)
        metric_low = aling_result.RMSDAfterRefinement

        if metric_low < 0 or low_actual_value < 0:
          best_iterations_num += (best_iterations_num - i)
          raise Exception

      except:
        metric_low = float("inf")
        substraction = 0

      try:
        if metric_high == float("inf"):
          raise Exception
        aling_result = align_pdb_file_1_in_2(pdb_simulated_high, original)
        print("Print 3", high_actual_value)
        print(aling_result.RMSDAfterRefinement)
        print(aling_result.PercentageAtomsAlignedAfterRefinement)
        print(aling_result.RawAlignmentScore)
        metric_high = aling_result.RMSDAfterRefinement

        if metric_high < 0:
          best_iterations_num += (best_iterations_num - i)
          raise Exception

      except:
        metric_high = float("inf")
        addition = 0

      if metric_low == float("inf") and metric_high == float("inf"):
        break
      elif metric_low <= metric_high and metric_low > 0:
        if metric_low <= last_metric:
          last_metric = metric_low
          level_to_return = low_actual_value
      else:
        if metric_high <= last_metric and metric_high > 0:
          last_metric = metric_high
          level_to_return = high_actual_value

      if metric_low <= break_sim or metric_high <= break_sim:
        break

      low_actual_value -= substraction
      high_actual_value += addition
      clean_file(low_simulate_path)
      clean_file(high_simulate_path)
      i += 1

    level = level_to_return

  free_file(low_simulate_path)
  free_file(high_simulate_path)
  free_file(simulate_path)

  return level


def get_mrc_to_pdb_aux(level, mrc_path, name_file):
  pdb_path = name_file
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


def mrc_to_pdb_test(mrc_path, pdb_result_path, original_path=None, threshold_mrc=0.0, clean=False):
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

  # execute_command("echo 1|../binaries/Situs_3.1/bin/map2map {0} {1}/tempPDB.situs".format(mrc_path, temp_dir))

  level = get_mrc_level(mrc_path)
  # Add flag, if original archive is here, find the best level, else, dont find best level

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
