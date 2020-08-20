'''
Created on Jul 19, 2020

@author: luis98
'''
import os
from subprocess import check_output, CalledProcessError
from tempfile import TemporaryFile
import shutil
from fit.fitMapResult import FitMapResult


def get_out(*args):
    with TemporaryFile() as t:
        try:
            out = check_output(args, stderr=t)
            return 0, out
        except CalledProcessError as e:
            t.seek(0)
            return e.returncode, t.read()


# Fit map0 into map1
def fit_map_in_map(map0_path, map1_path, path_exit_folder, attempts, map0_vector_move=None, map0_level=None, map1_level=None):
    path = "./temp_map"
    map0_exit_name = map0_path.split('/')[-1]
    map0_exit_name = map0_exit_name.split('.')[:-1]
    map0_exit_name = ".".join(map0_exit_name)
    map0_exit_name += "_fit.mrc"
    complete_exit_path = os.path.abspath(path_exit_folder)
    path_exit_folder = complete_exit_path + "/" + map0_exit_name

    if not os.path.exists(complete_exit_path):
        os.makedirs(complete_exit_path)

    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path)
    f = open(path + "/fit.cxc", "w+")

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

    f.write("fitmap #1 in_map #2 search {0} placement r\r\n".format(attempts))
    #f.write("fitmap #1 in_map #2\r\n")
    f.write("vop resample #1 ongrid #2 \r\n")
    f.write("save {0} #3\r\n".format(path_exit_folder))
    f.write("exit")
    f.close()

    map0_real_path = os.path.abspath(map0_path)
    map1_real_path = os.path.abspath(map1_path)
    commands_real_path = os.path.abspath(path + "/fit.cxc")

    _error, exit_binary_text = get_out("chimerax", "--nogui", map0_real_path, map1_real_path, commands_real_path)

    shutil.rmtree(path)
    text = exit_binary_text.decode("utf-8")
    # print(text)
    return FitMapResult(text, 'x')
