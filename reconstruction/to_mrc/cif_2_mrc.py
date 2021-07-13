import os

from general_utils.cif_utils import cif_to_pdb, get_chains_cif
from general_utils.temp_utils import gen_dir, free_dir
from general_utils.terminal_utils import get_out
from general_utils.pdb_utils import get_cube_pdb, move_pdb_center, get_atoms_of_list_pdb


def cif_to_mrc_chains(input_file, output_dir, chains=None, div_can=1):
  input_file = os.path.abspath(input_file)
  output_dir = os.path.abspath(output_dir)

  work_dir = gen_dir()
  pdb_path = os.path.join(work_dir, "pdbFile.pdb")
  cif_to_pdb(input_file, pdb_path)

  move_pdb_center(pdb_path)
  all_chains = get_chains_cif(input_file)
  information_atoms = get_atoms_of_list_pdb(pdb_path, all_chains)

  name_of_pdb = input_file.split('/')[-1]
  name_of_pdb = name_of_pdb.split('.')[:-1]
  name_of_pdb = ".".join(name_of_pdb)

  directory = output_dir + "/" + name_of_pdb
  if not os.path.exists(directory):
    os.makedirs(directory)

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

      for work_chain in chains[before_count:count]:
        lines += information_atoms[work_chain]

      final_text = "".join(lines)

      pdb_path = directory + "/" + name_of_pdb + "_" + ''.join(chains[before_count:count]) + ".pdb"
      f = open(pdb_path, "w+")
      f.write(final_text)
      f.close()

  free_dir(work_dir)
