import os
import math


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
    max_val = 2 * max(math.ceil(x_actual), y_actual, z_actual)
    max_val += 10
    return [max_val, max_val, max_val]
