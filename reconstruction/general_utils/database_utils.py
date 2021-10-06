import re
import shutil

import pymongo

import zlib, json, base64
import numpy as np
import os

from tqdm import tqdm
from Bio import SeqIO

from general_utils.cif_utils import get_chains_cif
from general_utils.download_utils import download_pdb, download_cif, download_pdb_fasta
from general_utils.mrc_uilts import mrc_to_pdb
from general_utils.pdb_utils import get_chains_pdb, make_pdb_dir_temp
from general_utils.temp_utils import gen_dir, free_dir, gen_file
from general_utils.workspace_utils import is_work_in_cluster
from miscellaneous_utils.database_mis import put_bigfile_db, get_dicc_chains_files_info, get_bigfile_db
from to_mrc.cif_2_mrc import cif_to_mrc_chains
from to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from networkx.readwrite import json_graph
from bson.json_util import dumps

database_name = "lcastillo_biostructures"
collection_name = "pdb_collection"
exists_mongo_db_var = None
valid_resolutions = [4, 6, 8, 10]
not_sequence = ["6v2y", "7lj6", "7lj8", "6v32", "6yz6", "7jml",
                "7lj8", "6v2y", "6v32", "7jml", "6yz6"]


# General mongo
def exists_mongo_db():
  # return False
  global exists_mongo_db_var
  if exists_mongo_db_var != None:
    return exists_mongo_db_var

  client = get_mongo_client()
  # Issue the serverStatus command and print the results
  try:
    client.server_info()
    exists_mongo_db_var = True
    return True
  except:
    exists_mongo_db_var = False
    return False


def get_mongo_client():
  if is_work_in_cluster():
    client = pymongo.MongoClient(host="11.0.0.21",
                                 port=27017,
                                 username='lcastilloAdmin',
                                 password='LPtYJpA3',
                                 authSource='lcastillo_biostructures')
  else:
    client = pymongo.MongoClient()
  return client


def clear_collection():
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    col.drop()


def save_collection(export_json_path):
  if os.path.exists(export_json_path):
    os.remove(export_json_path)

  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    cursor = col.find({}, no_cursor_timeout=True)
    with open(export_json_path, 'w') as file:
      file.write('[')

      flag_1 = True
      for document in cursor:
        if not flag_1:
          file.write(',')
        else:
          flag_1 = False
        file.write(dumps(document))
      file.write(']')


def load_collection(json_path):
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    with open(json_path) as f:
      file_data = json.load(f)
      if isinstance(file_data, list):
        with tqdm(total=len(file_data)) as pbar:
          for add_pdb in file_data:
            add_pdb.pop("_id")
            col.update_one(
              {"pdbID": add_pdb["pdbID"]},
              {"$setOnInsert": add_pdb},
              True
            )
            pbar.update(1)
      else:
        file_data.pop("_id")
        col.update_one(
          {"pdbID": file_data["pdbID"]},
          {"$setOnInsert": file_data},
          True
        )


def memory_use():
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    result = db.command("collstats", collection_name)
    final_size = result["totalIndexSize"] + result["storageSize"]
    return final_size / (1024 * 1024)
  return 0


# Creation
def create_get_pdb_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    if pdb_id in not_sequence:
      return False
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data is None:
      insert_pdb_information(db, col, pdb_id)
      pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)

    return pdb_data


def insert_pdb_information(db, col, pdb_id):
  global valid_resolutions
  from general_utils.pdb_utils import get_similar_pdb_struct, \
    get_similar_pdb_chain_structural, \
    get_similar_pdb_chain_sequential

  path_dir = gen_dir()
  is_pdb = True
  try:
    path_of_file = '{0}/{1}.pdb'.format(path_dir, pdb_id)
    download_pdb(pdb_id, path_of_file)
    chains = get_chains_pdb(path_of_file)
  except:
    is_pdb = False
    path_of_file = '{0}/{1}.cif'.format(path_dir, pdb_id)
    download_cif(pdb_id, path_of_file)
    chains = get_chains_cif(path_of_file)

  new_pdb = {}
  all_graphs = {}
  for resolution in valid_resolutions:
    key_mongo = pdb_id + "_" + str(resolution)
    dicc_graph = resolution_pdb_information(resolution, pdb_id, chains, path_dir, path_of_file, is_pdb)
    all_graphs[key_mongo] = dicc_graph

  new_pdb['all_graphs'] = put_bigfile_db(db, all_graphs)


  new_pdb['chains'] = put_bigfile_db(db, chains)
  all_sequences = get_online_sequences(pdb_id, chains)
  new_pdb['all_sequences'] = put_bigfile_db(db, all_sequences)

  map_pdb2cif = get_online_chain_map_pdb2cif(pdb_id, chains)
  new_pdb['map_pdb2cif'] = put_bigfile_db(db, map_pdb2cif)

  similars_pdbs = get_similar_pdb_struct(pdb_id, -1)
  new_pdb['similars_pdbs'] = put_bigfile_db(db, similars_pdbs)

  dicc_chain_struct = {}
  for chain in chains:
    dicc_chain_struct[chain] = get_similar_pdb_chain_structural(pdb_id, map_pdb2cif[chain], -1)
  new_pdb['similars_chain_struct'] = put_bigfile_db(db, dicc_chain_struct)

  dicc_chain_sequential = {}
  for chain in chains:
    dicc_chain_sequential[chain] = get_similar_pdb_chain_sequential(pdb_id, chain, -1, all_sequences[chain])
  new_pdb['similars_chain_sequential'] = put_bigfile_db(db, dicc_chain_sequential)

  chains_files_info = get_dicc_chains_files_info(pdb_id)
  new_pdb['chains_files_info'] = put_bigfile_db(db, chains_files_info)

  col.update_one(
    {"pdbID": pdb_id},
    {"$setOnInsert": new_pdb},
    True
  )

  free_dir(path_dir)
  return new_pdb


def resolution_pdb_information(resolution, pdb_id, chains, path_dir, path_of_file, is_pdb):
  graph, father_ZD = gen_graph_resolution_aux(chains, resolution, path_dir, pdb_id, path_of_file, is_pdb, True)

  for j in graph.nodes():
    graph.nodes[j]["zd_descriptors"] = graph.nodes[j]["zd_descriptors"].tolist()

  new_pdb_chain = {}

  new_pdb_chain["G"] = graph
  new_pdb_chain["OZD"] = father_ZD

  return new_pdb_chain


# Get online
# get chains online
def get_pdb_chains_online(pdb_id):
  path_dir = gen_dir()

  is_pdb = True
  try:
    path_of_file = '{0}/{1}.pdb'.format(path_dir, pdb_id)
    download_pdb(pdb_id, path_of_file)
    chains = get_chains_pdb(path_of_file)
  except:
    is_pdb = False
    path_of_file = '{0}/{1}.cif'.format(path_dir, pdb_id)
    download_cif(pdb_id, path_of_file)
    chains = get_chains_cif(path_of_file)

  free_dir(path_dir)
  return is_pdb, chains


def get_online_sequences(pdb, chains_check=None):
  if chains_check is None:
    chains_check = get_chains_pdb_db(pdb)

  if pdb in not_sequence:
    result = {}
    for i in chains_check:
      result[i] = ""
    return result

  work_dir = gen_dir()
  fasta_path = os.path.join(work_dir, "pdb.fasta")
  download_pdb_fasta(pdb, fasta_path, create_progress_bar=False)

  result = {}
  result_original = {}
  fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
  for fasta in fasta_sequences:
    sequence = str(fasta.seq)
    description = str(fasta.description)
    chains_sec = description.split("|")[1]

    chains_sec = chains_sec.replace("auths ", "")
    chains_sec = chains_sec.replace("auth ", "")

    chains_sec = chains_sec.replace("Chains ", "")
    chains_sec = chains_sec.replace("Chain ", "")

    chains_sec = chains_sec.replace(" ", "")

    # chains_sec = re.sub(r'\[[^)]*\]', '', chains_sec)

    # chains = chains_sec.split(",")
    # chains = re.split(r'\[*\]', chains_sec)
    chains = re.findall(r'[A-Za-z0-9]+\[[^]]*\]|[A-Za-z0-9]+', chains_sec)

    for chain in chains:
      list_temp = chain.split("[")
      chain_cif = list_temp[0]
      chain_cif = chain_cif.replace(",", "")

      result[chain_cif] = sequence

      if len(list_temp) > 1:
        chains_quivalent = list_temp[1].replace("],", "")
        chains_quivalent = chains_quivalent.replace("]", "")
        chains_quivalent = chains_quivalent.split(",")
        for i in chains_quivalent:
          result_original[i] = sequence
      else:
        result_original[chain_cif] = sequence

  free_dir(work_dir)

  origin_in = True
  for i in chains_check:
    if i not in result_original.keys():
      origin_in = False
      break

  if origin_in:
    return result_original

  cif_in = True
  for i in chains_check:
    if i not in result.keys():
      cif_in = False
      break

  if cif_in:
    return result

  raise ValueError("Cant not mapping sequence")


def get_online_chain_map_pdb2cif(pdb, chains):
  if pdb in not_sequence or not havePDBFile(pdb):
    result = {}
    for i in chains:
      result[i] = i
    return result

  work_dir = gen_dir()
  fasta_path = os.path.join(work_dir, "pdb.fasta")
  download_pdb_fasta(pdb, fasta_path, create_progress_bar=False)

  temp_result = {}
  fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
  for fasta in fasta_sequences:
    description = str(fasta.description)
    chains_sec = description.split("|")[1]

    chains_sec = chains_sec.replace("auths ", "")
    chains_sec = chains_sec.replace("auth ", "")

    chains_sec = chains_sec.replace("Chains ", "")
    chains_sec = chains_sec.replace("Chain ", "")

    chains_sec = chains_sec.replace(" ", "")

    chains = re.findall(r'[A-Za-z0-9]+\[[^]]*\]|[A-Za-z0-9]+', chains_sec)

    for chain in chains:
      list_temp = chain.split("[")
      chain_cif = list_temp[0]
      chain_cif = chain_cif.replace(",", "")

      temp_result[chain_cif] = []

      if len(list_temp) > 1:
        chains_quivalent = list_temp[1].replace("],", "")
        chains_quivalent = chains_quivalent.replace("]", "")
        chains_quivalent = chains_quivalent.split(",")
        for i in chains_quivalent:
          temp_result[chain_cif].append(i)
      else:
        temp_result[chain_cif].append(chain_cif)

  final_result = {}
  for i in temp_result.keys():
    for j in temp_result[i]:
      final_result[j] = i

  return final_result


def havePDBFile(pdb_id):
  work_local_dir = gen_dir()
  is_pdb = True
  try:
    path_of_file = '{0}/{1}.pdb'.format(work_local_dir, pdb_id)
    download_pdb(pdb_id, path_of_file)
    chains = get_chains_pdb(path_of_file)
  except:
    is_pdb = False
    path_of_file = '{0}/{1}.cif'.format(work_local_dir, pdb_id)
    download_cif(pdb_id, path_of_file)
    chains = get_chains_cif(path_of_file)

  free_dir(work_local_dir)
  return is_pdb


# Delete
def delete_pdb_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    _result = col.delete_one({'pdbID': pdb_id})


# General PDBS funtions
def get_all_archive_pdb():
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdbs_data = col.find(no_cursor_timeout=True)
    if pdbs_data.count() == 0:
      return []
    result = []
    for i in pdbs_data:
      result.append(i['pdbID'])
    return result

  else:
    return []


# General gets
def get_graph_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    all_graphs = get_bigfile_db(db, pdb_data['all_graphs'])
    key = pdb_id + "_" + resolution
    graph_result = all_graphs[key]["G"]
    return graph_result


def get_zd_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    all_graphs = get_bigfile_db(db, pdb_data['all_graphs'])
    key = pdb_id + "_" + resolution
    zd_complete = all_graphs[key]["OZD"]
    return zd_complete


def get_zd_chains_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    all_graphs = get_bigfile_db(db, pdb_data['all_graphs'])
    key = pdb_id + "_" + resolution
    graph = all_graphs[key]["G"]
    final_result = []

    for j in graph.nodes():
      final_result.append([j, np.array(graph.nodes[j]["zd_descriptors"])])

    return final_result


def get_zd_chain_pdb_db(pdb_id, chain, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    all_graphs = get_bigfile_db(db, pdb_data['all_graphs'])
    key = pdb_id + "_" + resolution
    graph = all_graphs[key]["G"]

    result = np.array(graph.nodes[chain]["zd_descriptors"])
    return result


def get_chains_pdb_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    chains_result = get_bigfile_db(db, pdb_data['chains'])

    return chains_result


def get_sequence_pdb_db(pdb_id, chain):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    all_sequences = get_bigfile_db(db, pdb_data['all_sequences'])
    result = all_sequences[chain]

    return result


def get_pdb2cif_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    pdb_data = create_get_pdb_db(pdb_id)
    client = get_mongo_client()
    db = client[database_name]

    map_pdb2cif = get_bigfile_db(db, pdb_data['map_pdb2cif'])
    return map_pdb2cif


# Aux funtions
def gen_graph_resolution_aux(chains, resolution, path_dir, pdb_id, path_file, is_pdb, return_OZD=False):
  if is_pdb:
    pdb_to_mrc_chains(return_OZD, False, resolution, path_file, path_dir, chains, len(chains))
  else:
    cif_to_mrc_chains(return_OZD, False, resolution, path_file, path_dir, chains, len(chains))

  local_path = path_dir + "/" + pdb_id
  segments = []
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_id, chain))

    segment = segments_graph_simulate[0]
    segment.id_segment = chain
    segment.textSimulatePDB = get_text_simulate_PDB('{0}/{1}_{2}.mrc'.format(local_path, pdb_id, chain))
    # segment.id_segment = chain
    segments.append(segment)
  graph = generate_graph(segments, 100, 0, 6, 1)
  if return_OZD:
    father_ZD = get_mrc_one('{0}/{1}.mrc'.format(local_path, pdb_id))[1].zd_descriptors
    free_dir(local_path)
    return graph, father_ZD
  else:
    free_dir(local_path)
    return graph

def get_text_simulate_PDB(mrc_file):
  file_tem = gen_file()
  mrc_to_pdb(mrc_file, file_tem, clean=True)

  result = None
  with open(file_tem) as f:
    lines = f.readlines()
    result = "".join(lines)
  os.remove(file_tem)
  return result


# client = get_mongo_client()
# Issue the serverStatus command and print the results
# serverStatusResult = db.command("serverStatus")
# pprint(serverStatusResult)
# print(exists_mongo_db())
# print(exists_mongo_db())
# graph = get_graph_pdb_db('100d', 4)
# result = get_zd_chains_pdb_db('100d', '2')
# result2 = get_zd_pdb_db('100d', 4)
# print(result2)
#
#
# print(get_chains_pdb_db('100d'))
# print(get_all_archive_pdb())
