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
from general_utils.pdb_utils import get_chains_pdb
from general_utils.temp_utils import gen_dir, free_dir
from general_utils.workspace_utils import is_work_in_cluster
from to_mrc.cif_2_mrc import cif_to_mrc_chains
from to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_graph.process_graph_utils import generate_graph
from process_mrc.generate import get_mrc_one
from networkx.readwrite import json_graph
from bson.json_util import dumps

ZIPJSON_KEY = 'base64(zip(o))'
database_name = "lcastillo_biostructures"
collection_name = "pdb_collection"
exists_mongo_db_var = None
valid_resolutions = [4, 6, 8, 10]
not_sequence = ["6v2y", "7lj6", "7lj8", "6v32", "6yz6", "7jml",
                "7lj8", "6v2y", "6v32", "7jml", "6yz6"]


def json_zip(j):
  j = {
    ZIPJSON_KEY: base64.b64encode(
      zlib.compress(
        json.dumps(j).encode('utf-8')
      )
    ).decode('ascii')
  }

  return j


def json_unzip(j, insist=True):
  try:
    assert (j[ZIPJSON_KEY])
    assert (set(j.keys()) == {ZIPJSON_KEY})
  except:
    if insist:
      raise RuntimeError("JSON not in the expected format {" + str(ZIPJSON_KEY) + ": zipstring}")
    else:
      return j

  try:
    j = zlib.decompress(base64.b64decode(j[ZIPJSON_KEY]))
  except:
    raise RuntimeError("Could not decode/unzip the contents")

  try:
    j = json.loads(j)
  except:
    raise RuntimeError("Could interpret the unzipped contents")

  return j


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


def get_graph_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      result_zip = pdb_data[str(resolution) + "A_G"]
      js_graph = json_unzip(result_zip)
      graph = json_graph.node_link_graph(js_graph)
      for j in graph.nodes():
        graph.nodes[j]["zd_descriptors"] = np.array(graph.nodes[j]["zd_descriptors"])
      return graph

    else:
      insert_pdb_information(col, pdb_id)
      return get_graph_pdb_db(pdb_id, resolution)

  else:
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

    graph = gen_graph_resolution_aux(chains, resolution, path_dir, pdb_id, path_of_file, is_pdb)
    free_dir(path_dir)
    return graph


def get_zd_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      result_zip = pdb_data[str(resolution) + "A_OZD"]
      js_zd = json_unzip(result_zip)
      list_js = json.loads(js_zd)
      final_result = np.array(list_js)
      return final_result

    else:
      insert_pdb_information(col, pdb_id)
      return get_zd_pdb_db(pdb_id, resolution)

  else:
    result = gen_zd_resolution_aux(resolution, pdb_id)
    return result


def get_zd_chains_pdb_db(pdb_id, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      result_zip = pdb_data[str(resolution) + "A_G"]
      js_graph = json_unzip(result_zip)
      graph = json_graph.node_link_graph(js_graph)
      final_result = []
      for j in graph.nodes():
        final_result.append([pdb_data["chains"][j - 1], np.array(graph.nodes[j]["zd_descriptors"])])
      return final_result

    else:
      insert_pdb_information(col, pdb_id)
      return get_zd_chains_pdb_db(pdb_id, resolution)

  else:
    result = gen_zd_chain_resolution_aux(resolution, pdb_id)
    return result


def get_chain_to_number_chain(pdb_id, chain):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      return pdb_data["chains"].index(chain) + 1


def get_zd_chain_pdb_db(pdb_id, chain, resolution):
  if resolution not in valid_resolutions:
    raise TypeError("Resolution is not valid, you can use: {0}".format(valid_resolutions))

  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      result_zip = pdb_data[str(resolution) + "A_G"]
      js_graph = json_unzip(result_zip)
      graph = json_graph.node_link_graph(js_graph)
      final_result = []
      for j in graph.nodes():
        final_result.append([pdb_data["chains"][j - 1], np.array(graph.nodes[j]["zd_descriptors"])])
      for k in final_result:
        if k[0] == chain:
          return k[1]

    else:
      insert_pdb_information(col, pdb_id)
      return get_zd_chains_pdb_db(pdb_id, resolution)

  else:
    result = gen_zd_chain_resolution_one_chain_aux(resolution, pdb_id, chain)
    return result


def get_chains_pdb_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      return pdb_data["chains"]

    else:
      insert_pdb_information(col, pdb_id)
      return get_chains_pdb_db(pdb_id)

  else:
    path_dir = gen_dir()
    download_pdb(pdb_id, '{0}/{1}.pdb'.format(path_dir, pdb_id))
    chains = get_chains_pdb('{0}/{1}.pdb'.format(path_dir, pdb_id))
    free_dir(path_dir)
    return chains


def get_sequence_pdb_db(pdb_id, chain):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      if 'all_sequences' in pdb_data.keys():
        dicc = pdb_data["all_sequences"]
        try:
          return dicc[chain]
        except:
          all_sequences = insert_sequences_db(pdb_id)
          return all_sequences[chain]
      else:
        insert_sequences_db(pdb_id)
        return get_sequence_pdb_db(pdb_id, chain)

    else:
      insert_pdb_information(col, pdb_id)
      return get_sequence_pdb_db(pdb_id, chain)

  else:
    all_sequences = get_online_sequences(pdb_id)
    return all_sequences[chain]


def insert_sequences_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    all_sequences = get_online_sequences(pdb_id)
    col.update_one(
      {"pdbID": pdb_id},
      {"$set": {'all_sequences': all_sequences}},
    )

    return all_sequences


def make_pdb_dir_temp(work_dir, pdb_id):
  complete_path = os.path.abspath(work_dir)
  dirs = os.listdir(complete_path)

  if pdb_id in dirs:
    if len(os.listdir(os.path.join(complete_path, pdb_id))) > 0:
      return
    else:
      os.rmdir(os.path.join(complete_path, pdb_id))

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

  if is_pdb:
    from to_mrc.pdb_2_mrc import pdb_to_chains_file
    pdb_to_chains_file(path_of_file, complete_path, chains, len(chains))
  else:
    from general_utils.cif_utils import cif_to_chains_pdb_files
    cif_to_chains_pdb_files(path_of_file, complete_path, chains, len(chains))


def get_db_chains_files_db(pdb_id, path):
  pdb_id = pdb_id.lower()
  make_pdb_dir_temp(path, pdb_id)
  # if exists_mongo_db():
  #   client = get_mongo_client()
  #   db = client[database_name]
  #   col = db[collection_name]
  #   pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
  #   if pdb_data != None:
  #     if 'chains_files_info' in pdb_data.keys():
  #       dicc = pdb_data["chains_files_info"]
  #
  #       path_result = os.path.join(path, pdb_id)
  #       shutil.rmtree(path_result, ignore_errors=True)
  #       os.mkdir(path_result)
  #
  #       for i in dicc.keys():
  #         f = open(os.path.join(path_result, i), 'w')
  #         f.write("".join(dicc[i]))
  #         f.close()
  #     else:
  #       dicc_add = get_dicc_chains_files_info(path, pdb_id)
  #
  #       try:
  #         col.update_one(
  #           {"pdbID": pdb_id},
  #           {"$set": {'chains_files_info': dicc_add}},
  #         )
  #       except:
  #         pass
  #   else:
  #     make_pdb_dir_temp(path, pdb_id)
  #     # insert_pdb_information(col, pdb_id)
  #     # return get_db_chains_files_db(pdb_id, path)
  #
  # else:
  #   make_pdb_dir_temp(path, pdb_id)


def get_dicc_chains_files_info(path, pdb_id):
  make_pdb_dir_temp(path, pdb_id)
  dicc_add = {}
  path_pdb = os.path.join(path, pdb_id)
  for chain_file in os.listdir(path_pdb):
    path_chain = os.path.join(path_pdb, chain_file)
    with open(path_chain) as f:
      dicc_add[os.path.basename(chain_file)] = f.readlines()
  return dicc_add


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


def get_pdb2cif_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    pdb_data = col.find_one({'pdbID': pdb_id}, no_cursor_timeout=True)
    if pdb_data != None:
      if 'map_pdb2cif' in pdb_data.keys():
        dicc = pdb_data["map_pdb2cif"]
        if dicc == {}:
          dicc = insert_map_pdb2cif_db(pdb_id, pdb_data["chains"])
        return dicc
      else:
        insert_map_pdb2cif_db(pdb_id, pdb_data["chains"])
        return get_pdb2cif_db(pdb_id)

    else:
      insert_pdb_information(col, pdb_id)
      return get_pdb2cif_db(pdb_id)

  else:
    dicc = get_online_chain_map_pdb2cif(pdb_id)
    return dicc


def insert_map_pdb2cif_db(pdb_id, chains):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    map_pdb2cif = get_online_chain_map_pdb2cif(pdb_id, chains)
    col.update_one(
      {"pdbID": pdb_id},
      {"$set": {'map_pdb2cif': map_pdb2cif}},
    )

    return map_pdb2cif


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


def delete_pdb_db(pdb_id):
  pdb_id = pdb_id.lower()
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    result = col.delete_one({'pdbID': pdb_id})


def insert_pdb_information(col, pdb_id):
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
  # new_pdb['pdbID'] = pdb_id
  new_pdb['chains'] = chains
  new_pdb['all_sequences'] = get_online_sequences(pdb_id, chains)
  new_pdb['map_pdb2cif'] = get_online_chain_map_pdb2cif(pdb_id, chains)

  for i in [4, 6, 8, 10]:
    graph, father_ZD = gen_graph_resolution_aux(chains, i, path_dir, pdb_id, path_of_file, is_pdb, True)
    for j in graph.nodes():
      graph.nodes[j]["zd_descriptors"] = graph.nodes[j]["zd_descriptors"].tolist()
    g_json = json_graph.node_link_data(graph)
    new_pdb[str(i) + "A_G"] = json_zip(g_json)
    new_pdb[str(i) + "A_OZD"] = json_zip(json.dumps(father_ZD.tolist()))

  col.update_one(
    {"pdbID": pdb_id},
    {"$setOnInsert": new_pdb},
    True
  )

  free_dir(path_dir)
  return new_pdb


def gen_graph_resolution_aux(chains, resolution, path_dir, pdb_id, path_file, is_pdb, return_OZD=False):
  if is_pdb:
    pdb_to_mrc_chains(return_OZD, False, resolution, path_file, path_dir, chains, len(chains))
  else:
    cif_to_mrc_chains(return_OZD, False, resolution, path_file, path_dir, chains, len(chains))

  local_path = path_dir + "/" + pdb_id
  con_id_segment = 1
  segments = []
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_id, chain))
    segment = segments_graph_simulate[0]
    segment.id_segment = con_id_segment
    # segment.id_segment = chain
    segments.append(segment)
    con_id_segment += 1
  graph = generate_graph(segments, 100, 0, 6, 1)
  if return_OZD:
    father_ZD = get_mrc_one('{0}/{1}.mrc'.format(local_path, pdb_id))[1].zd_descriptors
    free_dir(local_path)
    return graph, father_ZD
  else:
    free_dir(local_path)
    return graph


def gen_zd_chain_resolution_aux(resolution, pdb_id):
  path_dir = gen_dir()

  try:
    path_of_pdb = '{0}/{1}.pdb'.format(path_dir, pdb_id)
    download_pdb(pdb_id, path_of_pdb)
    chains = get_chains_pdb(path_of_pdb)
    pdb_to_mrc_chains(False, False, resolution, path_of_pdb, path_dir, chains, len(chains))
  except:
    path_of_cif = '{0}/{1}.cif'.format(path_dir, pdb_id)
    download_cif(pdb_id, path_of_cif)
    chains = get_chains_cif(path_of_cif)
    cif_to_mrc_chains(False, False, resolution, path_of_cif, path_dir, chains, len(chains))

  local_path = path_dir + "/" + pdb_id
  result = []
  for chain in chains:
    segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_id, chain))
    segment = segments_graph_simulate[0]
    result.append([chain, segment.zd_descriptors])

  free_dir(path_dir)
  return result


def gen_zd_chain_resolution_one_chain_aux(resolution, pdb_id, chain):
  path_dir = gen_dir()

  try:
    path_of_pdb = '{0}/{1}.pdb'.format(path_dir, pdb_id)
    download_pdb(pdb_id, path_of_pdb)
    # chains = get_chains_pdb(path_of_pdb)
    chains = [chain]
    pdb_to_mrc_chains(False, False, resolution, path_of_pdb, path_dir, chains, len(chains))
  except:
    path_of_cif = '{0}/{1}.cif'.format(path_dir, pdb_id)
    download_cif(pdb_id, path_of_cif)
    # chains = get_chains_cif(path_of_cif)
    chains = [chain]
    cif_to_mrc_chains(False, False, resolution, path_of_cif, path_dir, chains, len(chains))

  local_path = path_dir + "/" + pdb_id
  result = []

  segments_graph_simulate, _ = get_mrc_one('{0}/{1}_{2}.mrc'.format(local_path, pdb_id, chain))
  segment = segments_graph_simulate[0]
  result.append([chain, segment.zd_descriptors])

  free_dir(path_dir)
  return result


def gen_zd_resolution_aux(resolution, pdb_id):
  path_dir = gen_dir()

  try:
    path_of_pdb = '{0}/{1}.pdb'.format(path_dir, pdb_id)
    download_pdb(pdb_id, path_of_pdb)
    pdb_to_mrc_chains(True, False, resolution, path_of_pdb, path_dir)
  except:
    path_of_cif = '{0}/{1}.cif'.format(path_dir, pdb_id)
    download_cif(pdb_id, path_of_cif)
    cif_to_mrc_chains(True, False, resolution, path_of_cif, path_dir)

  local_path = path_dir + "/" + pdb_id
  father_ZD = get_mrc_one('{0}/{1}.mrc'.format(local_path, pdb_id))[1].zd_descriptors
  free_dir(path_dir)
  return father_ZD


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


def get_dicc_pdbs_can_chains():
  result = {}
  if exists_mongo_db():
    client = get_mongo_client()
    db = client[database_name]
    col = db[collection_name]
    cursor = col.find({}, no_cursor_timeout=True)
    for pdb_data in cursor:
      chains_len = len(pdb_data["chains"])
      if not (chains_len in result.keys()):
        result[chains_len] = [pdb_data["pdbID"]]
      else:
        result[chains_len].append(pdb_data["pdbID"])

  return result


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
