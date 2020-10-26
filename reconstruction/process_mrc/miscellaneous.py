import shutil
import numpy as np
import os
from general_utils.terminal_utils import get_out
from general_utils.string_utils import get_float_between_ss, get_float_value


def get_center_point(index_segments, segments, filter_val):
    points_pos = None

    for i in index_segments:
        for j in segments:
            if j.id_segment == i:
                actual_segment = j
                temp_points = np.where(actual_segment.mask > filter_val)
                if points_pos is None:
                    points_pos = np.column_stack(temp_points)
                else:
                    points_pos = np.concatenate((points_pos, np.column_stack(temp_points)), axis=0)

    center_point = np.mean(points_pos, axis=0)
    return [int(center_point[0]), int(center_point[1]), int(center_point[2])]


def get_cube_len_angstrom(map_path):
    map_real_path = os.path.abspath(map_path)
    path = "./temp_map_center"
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

    shutil.rmtree(path)

    return [x, y, z]


def get_mass_angstrom(map_path):
    map_real_path = os.path.abspath(map_path)
    path = "./temp_map_mass"
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
        raise Exception("Error on try to get mass")
    text = exit_binary_text
    mass = get_float_between_ss(text, "Enclosed volume for surface (#1.1) =", "\n")
    shutil.rmtree(path)
    return mass


def get_mrc_level(map_path):
    map_real_path = os.path.abspath(map_path)
    path = "./temp_map_mass"
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
        raise Exception("Error on try to get mass")
    text = exit_binary_text
    level = get_float_between_ss(text, "at level", ",")
    shutil.rmtree(path)
    return level
