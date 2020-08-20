import os
import math
import shutil
import urllib
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


def download_pdb(id_code, exit_path):
    exit_path_full = os.path.abspath(exit_path)
    file_url = "https://files.rcsb.org/download/{0}.pdb".format(str(id_code))
    urllib.request.urlretrieve(file_url, exit_path_full)


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
    max_val = max(math.ceil(x_actual), y_actual, z_actual)
    max_val += 10
    return [max_val, max_val, max_val]


def pdb_to_mrc_chains(create_original, verbose, resolution, input_file, output_dir, chains=None, div_can=1,
                      cube_dimentions=None):

    if cube_dimentions is None:
        cube_dimentions = get_cube_pdb(input_file)

    input_file = os.path.abspath(input_file)
    output_dir = os.path.abspath(output_dir)

    name_of_pdb = input_file.split('/')[-1]
    name_of_pdb = name_of_pdb.split('.')[:-1]
    name_of_pdb = ".".join(name_of_pdb)
    out_file = "{0}.mrc".format(name_of_pdb)

    directory = output_dir + "/" + name_of_pdb
    if not os.path.exists(directory):
        os.makedirs(directory)
    complete_file_path = directory + "/" + str(out_file)

    if create_original:
        _exit, output = get_out('e2pdb2mrc.py', '-R', str(resolution), '-B',
                                '{0},{1},{2}'.format(cube_dimentions[0], cube_dimentions[1], cube_dimentions[2]),
                                str(input_file), complete_file_path)

        if verbose:
            print(output)

    if chains is not None:
        div_can = min(div_can, len(chains))
        count = 0
        before_count = 0
        increase = len(chains) // div_can

        while count < div_can:
            before_count = count
            if count + increase > len(chains):
                count = len(chains)
            else:
                count += increase

            lines = []

            with open(input_file) as origin_file:
                actual_chain = ''
                for line in origin_file:
                    if line[0:4] == "ATOM":
                        if actual_chain == '':
                            actual_chain = line[21:22]
                        elif actual_chain != line[21:22]:
                            actual_chain = line[21:22]

                        if actual_chain in chains[before_count:count]:
                            lines.append(line)

                final_text = "".join(lines)

                pdb_path = directory + "/" + name_of_pdb + "_" + ''.join(chains[before_count:count]) + ".pdb"
                f = open(pdb_path, "w+")
                f.write(final_text)
                f.close()

                exit_mrc_path = directory + "/" + name_of_pdb + "_" + ''.join(chains[before_count:count]) + ".mrc"

                _exit, output = get_out('e2pdb2mrc.py', '-R', '-B', '{0},{1},{2}'.format(cube_dimentions[0],
                                                                                         cube_dimentions[1],
                                                                                         cube_dimentions[2]),
                                        str(resolution), str(pdb_path), exit_mrc_path)

                os.remove(pdb_path)

                if verbose:
                    print(output.decode("utf-8"))
