import shutil

import numpy as np
import os
from subprocess import check_output, CalledProcessError
from tempfile import TemporaryFile


def get_out(*args):
    with TemporaryFile() as t:
        try:
            out = check_output(args, stderr=t)
            return 0, out
        except CalledProcessError as e:
            t.seek(0)
            return e.returncode, t.read()


def get_center_point(index_segments, segments, filter):
    x_pos = np.array([])
    y_pos = np.array([])
    z_pos = np.array([])

    con = 0
    for i in index_segments:
        for j in segments:
            if j.id_segment == i:
                actual_segment = j
                temp_points = np.where(actual_segment.mask > filter)
                x_pos = np.concatenate((x_pos, temp_points[0]), axis=0)
                y_pos = np.concatenate((y_pos, temp_points[1]), axis=0)
                z_pos = np.concatenate((z_pos, temp_points[2]), axis=0)
                break
    return [int(np.median(x_pos)), int(np.median(y_pos)), int(np.median(z_pos))]


def transform_points_sscale_aux(point, shape, shape_other, pos):
    if max(shape[pos], shape_other[pos]) == shape[pos]:
        return point[pos]
    else:
        return point[pos] + ((shape_other[pos] - shape[pos]) // 2)


def transform_points_sscale(point1, point2, shape1, shape2):
    newp1 = [0, 0, 0]
    newp2 = [0, 0, 0]

    newp1[0] = transform_points_sscale_aux(point1, shape1, shape2, 0)
    newp1[1] = transform_points_sscale_aux(point1, shape1, shape2, 1)
    newp1[2] = transform_points_sscale_aux(point1, shape1, shape2, 2)

    newp2[0] = transform_points_sscale_aux(point2, shape2, shape1, 0)
    newp2[1] = transform_points_sscale_aux(point2, shape2, shape1, 1)
    newp2[2] = transform_points_sscale_aux(point2, shape2, shape1, 2)

    return newp1, newp2


def get_vector_move_1to2(point1, point2):
    result = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
    return result


def get_float_betwen_ss(tex, a, b):
    pos_a = tex.find(a)
    tex = tex[pos_a + len(a):-1]
    pos_b = tex.find(b)

    return float(tex[0:pos_b])


def get_float_value(tex, a, b):
    r = 0
    r = get_float_betwen_ss(tex, a, b)
    return r


def get_cube_len(map_path):
    map_real_path = os.path.abspath(map_path)
    path = "./temp_map_center"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path)
    f = open(path + "/fit.cxc", "w+")
    f.write("volume #1 dumpHeader true \r\n")
    f.write("exit")
    f.close()

    commands_real_path = os.path.abspath(path + "/fit.cxc")

    _error, exit_binary_text = get_out("chimerax", "--nogui", map_real_path, commands_real_path)
    text = exit_binary_text.decode("utf-8")
    x = get_float_value(text, 'xlen =', '\n')
    y = get_float_value(text, 'ylen =', '\n')
    z = get_float_value(text, 'zlen =', '\n')

    shutil.rmtree(path)

    return [x, y, z]


def chance_basec(point, actual_shape, new_shape):
    result = [0, 0, 0]

    result[0] = new_shape[0] * (point[0] / actual_shape[0])
    result[1] = new_shape[1] * (point[1] / actual_shape[1])
    result[2] = new_shape[2] * (point[2] / actual_shape[2])

    return result
