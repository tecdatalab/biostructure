'''
Created on Jul 19, 2020

@author: luis98
'''
import os
import shutil

from fit.fit_result_chimeraX import FitMapResult
from general_utils.terminal_utils import get_out


def create_execute_file(path, map0_level, map1_level, map0_vector_move, map1_vector_move, attempts, path_exit_folder,
                        is_global):
    f = open(path + "/fit.cxc", "w+")
    f.write("volume #1 origin 0,0,0 \r\n")
    f.write("volume #2 origin 0,0,0 \r\n")

    if map0_level is not None:
        # f.write("volume #1 level " + str(map0_level).replace(".", ",") + "\r\n")
        f.write("volume #1 level " + str(map0_level) + "\r\n")

    if map1_level is not None:
        # f.write("volume #2 level " + str(map1_level).replace(".", ",") + "\r\n")
        f.write("volume #2 level " + str(map1_level) + "\r\n")

    if map0_vector_move is not None:
        f.write("move x {0} models #1".format(map0_vector_move[0]) + "\r\n")
        f.write("move y {0} models #1".format(map0_vector_move[1]) + "\r\n")
        f.write("move z {0} models #1".format(map0_vector_move[2]) + "\r\n")

    if map1_vector_move is not None:
        f.write("move x {0} models #2".format(map1_vector_move[0]) + "\r\n")
        f.write("move y {0} models #2".format(map1_vector_move[1]) + "\r\n")
        f.write("move z {0} models #2".format(map1_vector_move[2]) + "\r\n")
    if is_global:
        # f.write("fitmap #1 in_map #2 search {0} placement r\r\n".format(attempts))
        f.write("fitmap #1 in_map #2 metric correlation search {0} placement r\r\n".format(attempts))
    else:
        # f.write("fitmap #1 in_map #2 metric correlation \r\n")
        f.write("fitmap #1 in_map #2 metric correlation \r\n")
    # f.write("fitmap #1 in_map #2\r\n")
    f.write("vop resample #1 ongrid #2 \r\n")
    f.write("save {0} #3\r\n".format(path_exit_folder))
    f.write("volume mask #1 surface #2\r\n")
    f.write("measure volume #4\r\n")
    f.write("exit")
    f.close()


# Fit map0 into map1
def fit_map_in_map(map0_path, map1_path, path_exit_folder, attempts, map0_vector_move=None, map1_vector_move=None,
                   map0_level=None, map1_level=None):
    path = "./temp_map"
    map0_exit_name = map0_path.split('/')[-1]
    map0_exit_name = map0_exit_name.split('.')[:-1]
    map0_exit_name = ".".join(map0_exit_name)
    map0_exit_name += "_fit.mrc"
    complete_exit_path = os.path.abspath(path_exit_folder)
    path_exit_folder = complete_exit_path + "/" + map0_exit_name

    if not os.path.exists(complete_exit_path):
        os.makedirs(complete_exit_path)

    if os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)
    create_execute_file(path, map0_level, map1_level, map0_vector_move, map1_vector_move, attempts, path_exit_folder,
                        True)

    map0_real_path = os.path.abspath(map0_path)
    map1_real_path = os.path.abspath(map1_path)
    commands_real_path = os.path.abspath(path + "/fit.cxc")

    error, exit_binary_text = get_out("chimerax", "--nogui", map0_real_path, map1_real_path, commands_real_path)
    if error == 0:
        text = exit_binary_text
        if text.find('Found 0 fits') != -1:
            create_execute_file(path, map0_level, map1_level, map0_vector_move, map1_vector_move, attempts,
                                path_exit_folder,
                                False)
            error, exit_binary_text = get_out("chimerax", "--nogui", map0_real_path, map1_real_path, commands_real_path)
            text = exit_binary_text

        shutil.rmtree(path)
        # print(text)
        return FitMapResult(text, 'x')
    else:
        raise Exception("Error in fit map in map, of map1 {0} and map2 {1}, exit are: {2}".format(map0_real_path,
                                                                                                  map1_real_path,
                                                                                                  exit_binary_text.
                                                                                                  decode("utf-8")))
