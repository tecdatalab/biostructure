import os

from pymol import cmd

from general_utils.pdb_utils import move_pdb_center, get_atoms_of_list_pdb, get_pdb_chain_sequence, \
  get_chain_mmcif_2_pdb
from general_utils.temp_utils import gen_dir, free_dir
import numpy as np
import biotite.structure.io.pdbx as pdbx
import biotite


def get_chains_cif(input_file):
  input_file = os.path.abspath(input_file)

  work_dir = gen_dir()
  pdb_path = os.path.join(work_dir, "pdbFile.pdb")
  cif_to_pdb(input_file, pdb_path)
  ok_chains = get_chain_mmcif_2_pdb(pdb_path)

  # Read file
  file = pdbx.PDBxFile()
  file.read(input_file)
  # Get 'pdbx_struct_assembly_gen' category as dictionary
  assembly_dict = file["pdbx_struct_assembly_gen"]

  result = []
  if type(assembly_dict["asym_id_list"]) == str:
    result = assembly_dict["asym_id_list"].split(",")
  else:
    for asym_id_list in assembly_dict["asym_id_list"]:
      chain_ids = asym_id_list.split(",")
      for chain in chain_ids:
        if chain not in result:
          chain = chain.replace(" ", "")
          if not (chain == "" or \
                  chain is None or \
                  chain == " " or \
                  chain == '' or \
                  chain == ' '):
            result.append(chain)

  for chain in ok_chains:
    if chain not in result:
      error_msg = "Error chain:" + chain + ", error parse chains in cif, input_file: " + input_file + ", result: " + str(
        result) + ", ok_chains: " + str(ok_chains)

      # name_pdb_id = os.path.splitext(os.path.basename(input_file))[0]
      # from general_utils.database_utils import delete_pdb_db
      # delete_pdb_db(name_pdb_id)
      # from general_utils.database_utils import get_chains_pdb_db
      # get_chains_pdb_db(name_pdb_id)
      raise ValueError(error_msg)

  for chain in result[:]:
    if chain not in ok_chains:
      result.remove(chain)

  free_dir(work_dir)
  return result


def cif_to_pdb(mmCIF_path, exit_pdb_path):
  # Restar
  cmd.reinitialize()
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


def cif_to_chains_pdb_files(input_file, output_dir, chains=None, div_can=1):
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
