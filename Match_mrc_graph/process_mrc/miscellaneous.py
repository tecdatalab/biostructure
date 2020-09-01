import shutil
import numpy as np
import os
from general_utils.terminal_utils import get_out
from general_utils.string_utils import get_float_between_ss, get_float_value


def get_center_point(index_segments, segments, filter_val):
    x_pos = np.array([])
    y_pos = np.array([])
    z_pos = np.array([])

    for i in index_segments:
        for j in segments:
            if j.id_segment == i:
                actual_segment = j
                temp_points = np.where(actual_segment.mask > filter_val)
                x_pos = np.concatenate((x_pos, temp_points[0]), axis=0)
                y_pos = np.concatenate((y_pos, temp_points[1]), axis=0)
                z_pos = np.concatenate((z_pos, temp_points[2]), axis=0)

    return [int(np.median(x_pos)), int(np.median(y_pos)), int(np.median(z_pos))]


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
