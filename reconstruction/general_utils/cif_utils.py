import os

from pymol import cmd

from general_utils.pdb_utils import move_pdb_center, get_atoms_of_list_pdb, get_pdb_chain_sequence
from general_utils.temp_utils import gen_dir, free_dir


def get_chains_cif(input_file):
  input_file = os.path.abspath(input_file)
  list_result = []
  with open(input_file) as origin_file:
    actual_chain = ''
    for line in origin_file:
      if line[0:4] == "ATOM":
        if actual_chain == '':
          actual_chain = line[28:31].replace(" ", "")
        elif actual_chain != line[28:31].replace(" ", ""):
          if actual_chain not in list_result:
            list_result.append(actual_chain)
          actual_chain = line[28:31].replace(" ", "")

    if actual_chain not in list_result:
      list_result.append(actual_chain)
  return list_result


def cif_to_pdb(mmCIF_path, exit_pdb_path):
  exit_pdb_path = os.path.abspath(exit_pdb_path)
  input_file = os.path.abspath(mmCIF_path)
  pdb_name = os.path.basename(input_file).split('.')[0]
  cmd.load(input_file, pdb_name)
  cmd.save(exit_pdb_path)
  # cmd.quit()


def get_cif_chain_sequence(cif_path, pdb_name, chain):
  cif_path = os.path.abspath(cif_path)

  # from Bio import SeqIO
  #
  # for record in SeqIO.parse(cif_path, "cif-seqres"):
  #   print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
  #   print(record.dbxrefs)

  cif_path = os.path.abspath(cif_path)

  work_dir = gen_dir()
  pdb_path = os.path.join(work_dir, "pdbFile.pdb")
  cif_to_pdb(cif_path, pdb_path)

  move_pdb_center(pdb_path)
  all_chains = get_chains_cif(cif_path)
  information_atoms = get_atoms_of_list_pdb(pdb_path, all_chains)

  final_text = "".join(information_atoms[chain])

  pdb_chain_path = work_dir + "/" + pdb_name + "_Chain.pdb"
  f = open(pdb_chain_path, "w+")
  f.write(final_text)
  f.close()

  sequence = get_pdb_chain_sequence(pdb_chain_path, pdb_name, "A")

  free_dir(work_dir)
  return sequence
