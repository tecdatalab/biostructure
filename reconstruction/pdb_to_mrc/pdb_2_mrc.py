import os
from general_utils.terminal_utils import get_out
from pdb_to_mrc.miscellaneous import get_cube_pdb


def pdb_to_mrc_chains(create_original, verbose, resolution, input_file, output_dir, chains=None, div_can=1,
                      cube_dimensions=None):
    if cube_dimensions is None:
        cube_dimensions = get_cube_pdb(input_file)

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
                                '{0},{1},{2}'.format(cube_dimensions[0], cube_dimensions[1], cube_dimensions[2]),
                                str(input_file), complete_file_path)

        if verbose:
            print(output.decode("utf-8"))

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

                _exit, output = get_out('e2pdb2mrc.py', '-R', str(resolution), '-B',
                                        '{0},{1},{2}'.format(cube_dimensions[0], cube_dimensions[1],
                                                             cube_dimensions[2]), str(pdb_path), exit_mrc_path)

                os.remove(pdb_path)

                if verbose:
                    print(output.decode("utf-8"))
