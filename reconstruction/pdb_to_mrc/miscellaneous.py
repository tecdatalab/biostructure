import os
import math
import numpy as np

from general_utils.string_utils import change_string


def get_chains(input_file):
  input_file = os.path.abspath(input_file)
  list_result = []
  with open(input_file) as origin_file:
    actual_chain = ''
    for line in origin_file:
      if line[0:4] == "ATOM":
        if actual_chain == '':
          actual_chain = line[21:22]
        elif actual_chain != line[21:22]:
          list_result.append(actual_chain)
          actual_chain = line[21:22]
    list_result.append(actual_chain)
  list_result = list(dict.fromkeys(list_result))
  return list_result


def get_cube_pdb(input_file):
  input_file = os.path.abspath(input_file)
  x_actual = 0.0
  y_actual = 0.0
  z_actual = 0.0
  with open(input_file) as origin_file:
    for line in origin_file:
      if line[0:4] == "ATOM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        x_actual = max(x, x_actual)
        y_actual = max(y, y_actual)
        z_actual = max(z, z_actual)
  max_val = 2 * max(math.ceil(x_actual), math.ceil(y_actual), math.ceil(z_actual))
  max_val += 10
  return [max_val, max_val, max_val]


def move_pdb_center(pdb_path):
  input_file = os.path.abspath(pdb_path)
  atoms = []

  with open(input_file) as origin_file:
    for line in origin_file:
      if line[0:4] == "ATOM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atoms.append([x,y,z])

  center = np.mean(atoms, axis=0)

  all_lines = []
  with open(input_file) as origin_file2:
    for line in origin_file2:
      if line[0:4] == "ATOM":
        x = float(line[30:38]) - center[0]
        y = float(line[38:46]) - center[1]
        z = float(line[46:54]) - center[2]

        text_x = "{:8.3f}".format(x)
        text_y = "{:8.3f}".format(y)
        text_z = "{:8.3f}".format(z)

        line = change_string(30, 37, line, text_x)
        line = change_string(38, 45, line, text_y)
        line = change_string(46, 53, line, text_z)
      all_lines.append(line)

  final_text = "".join(all_lines)
  os.remove(pdb_path)
  exit_file = open(pdb_path, "w")
  exit_file.write(final_text)
  exit_file.close()

