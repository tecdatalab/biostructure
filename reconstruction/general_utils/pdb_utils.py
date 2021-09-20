import itertools
import json
import math
import pickle
import time
import random
import os
import pandas as pd
from ftplib import FTP

import numpy as np
import requests
import re

from general_utils.string_utils import change_string

from general_utils.temp_utils import gen_dir, free_dir
from models.PDBsAlingResult import PDBsAlingResult

adn_arn_online_list = []


def get_similar_pdb_struct(pdb_name, can=10):
  search_request = {
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {
          "type": "group",
          "logical_operator": "and",
          "nodes": [
            {
              "type": "group",
              "logical_operator": "and",
              "nodes": [
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "DNA-RNA"
                  }
                },
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "RNA"
                  }
                },
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "DNA"
                  }
                }
              ],
              "label": "input-group"
            },
            {
              "type": "terminal",
              "service": "text",
              "parameters": {
                "operator": "in",
                "negation": True,
                "value": list(map(lambda x: x.upper(), get_pdb_no_work())),
                "attribute": "rcsb_entry_container_identifiers.entry_id"
              }
            }
          ]
        },
        {
          "type": "terminal",
          "service": "structure",
          "parameters": {
            "value": {
              "entry_id": pdb_name.upper(),
              "assembly_id": "1"
            },
            "operator": "relaxed_shape_match"
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "return_all_hits": True,
      "scoring_strategy": "combined",
      "sort": [
        {
          "sort_by": "score",
          "direction": "desc"
        }
      ]
    }
  }
  json_dump = json.dumps(search_request)
  url_get = 'https://search.rcsb.org/rcsbsearch/v1/query?json={0}'.format(json_dump)

  status_code = 500
  while status_code != 200 and status_code != 204 and status_code != 400:
    response = requests.get(url_get)
    status_code = response.status_code
    time.sleep(random.randint(2, 15))
    if status_code != 200 and status_code != 204:
      if status_code == 500 and response.text.find("ERROR") != -1:
        return []
      print("Check why\n",
            "Only pdb",
            pdb_name,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)
    elif status_code == 204:
      if response.text == "":
        return []
      print("Check why\n",
            "Only pdb",
            pdb_name,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)

  if response.status_code == 204 or response.status_code == 400:
    return []
  var_result = json.loads(response.text)
  result = []
  # Make sure that only pdb appear that we can process
  chain_score = {}
  dirty_list = []
  for i in var_result["result_set"]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      dirty_list.append(i["identifier"].split("-")[0].lower())
      chain_score[i["identifier"].split("-")[0].lower()] = i["score"]

  clean_list = np.intersect1d(np.array(dirty_list), np.array(get_all_pdb_work())).tolist()
  for i in clean_list[:can]:
    result.append([i, chain_score[i]])

  return result


def get_similar_pdb_chain_structural(pdb_name, chain, can=10):
  from general_utils.database_utils import get_pdb2cif_db
  chain_searh = get_pdb2cif_db(pdb_name)[chain]
  search_request = {
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {
          "type": "group",
          "logical_operator": "and",
          "nodes": [
            {
              "type": "group",
              "logical_operator": "and",
              "nodes": [
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "DNA-RNA"
                  }
                },
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "RNA"
                  }
                },
                {
                  "type": "terminal",
                  "service": "text",
                  "parameters": {
                    "attribute": "struct_keywords.pdbx_keywords",
                    "negation": True,
                    "operator": "contains_phrase",
                    "value": "DNA"
                  }
                }
              ],
              "label": "input-group"
            },
            {
              "type": "terminal",
              "service": "text",
              "parameters": {
                "operator": "in",
                "negation": True,
                "value": list(map(lambda x: x.upper(), get_pdb_no_work())),
                "attribute": "rcsb_entry_container_identifiers.entry_id"
              }
            }
          ]
        },
        {
          "type": "terminal",
          "service": "structure",
          "parameters": {
            "value": {
              "entry_id": pdb_name.upper(),
              "asym_id": chain_searh
            },
            "operator": "relaxed_shape_match"
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "return_all_hits": True,
      "scoring_strategy": "combined",
      "sort": [
        {
          "sort_by": "score",
          "direction": "desc"
        }
      ]
    }
  }

  json_dump = json.dumps(search_request)
  url_get = 'https://search.rcsb.org/rcsbsearch/v1/query?json={0}'.format(json_dump)

  status_code = 500
  while status_code != 200 and status_code != 204 and status_code != 400:
    response = requests.get(url_get)
    status_code = response.status_code
    time.sleep(random.randint(2, 15))
    if status_code != 200 and status_code != 204:
      if status_code == 500 and response.text.find("ERROR") != -1:
        return []
      print("Check why\n",
            "Only chain struct",
            pdb_name,
            "\n",
            chain + "=>" + chain_searh,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)
    elif status_code == 204:
      if response.text == "":
        return []
      print("Check why\n",
            "Only chain struct",
            pdb_name,
            "\n",
            chain + "=>" + chain_searh,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)
  if response.status_code == 204 or response.status_code == 400:
    return []
  var_result = json.loads(response.text)
  result = []

  # Make sure that only pdb appear that we can process
  chain_score = {}
  dirty_list = []
  for i in var_result["result_set"]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      dirty_list.append(i["identifier"].split("-")[0].lower())
      chain_score[i["identifier"].split("-")[0].lower()] = i["score"]

  clean_list = np.intersect1d(np.array(dirty_list), np.array(get_all_pdb_work())).tolist()
  for i in clean_list[:can]:
    result.append([i, chain_score[i]])

  return result


def type_sequence(sequence):
  math_protein = '^[ARNDCEQGHILKMFPSTWYVX]+$'
  math_dna = '^[ACTG]+$'
  math_rna = '^[ACUG]+$'

  if re.match(math_dna, sequence):
    return "pdb_dna_sequence"
  elif re.match(math_rna, sequence):
    return "pdb_rna_sequence"
  elif re.match(math_protein, sequence):
    return "pdb_protein_sequence"
  else:
    raise NameError('Can not get type of sequence')


def get_similar_pdb_chain_sequential(pdb_name, chain, can=10):
  from general_utils.database_utils import get_sequence_pdb_db
  sequence = get_sequence_pdb_db(pdb_name, chain)

  try:
    type_s = type_sequence(sequence)
  except:
    return []

  if len(sequence) < 10 or sequence == len(sequence) * 'X':
    return []

  search_request = {
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {
          "type": "group",
          "logical_operator": "and",
          "nodes": [
            {
              "type": "terminal",
              "service": "text",
              "parameters": {
                "operator": "contains_phrase",
                "negation": True,
                "value": "DNA,RNA,DNA-RNA",
                "attribute": "struct_keywords.pdbx_keywords"
              }
            }
          ]
        },
        {
          "type": "terminal",
          "service": "sequence",
          "parameters": {
            "evalue_cutoff": 0.1,
            "identity_cutoff": 0,
            "target": type_s,
            "value": sequence
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "return_all_hits": True,
      "scoring_strategy": "combined",
      "sort": [
        {
          "sort_by": "score",
          "direction": "desc"
        }
      ]
    }
  }

  json_dump = json.dumps(search_request)
  url_get = 'https://search.rcsb.org/rcsbsearch/v1/query?json={0}'.format(json_dump)

  status_code = 500
  while status_code != 200 and status_code != 204 and status_code != 400:
    response = requests.get(url_get)
    status_code = response.status_code
    time.sleep(random.randint(2, 15))
    if status_code != 200 and status_code != 204:
      if status_code == 500 and response.text.find("ERROR") != -1:
        return []
      print("Check why\n",
            "Only chain sequence",
            pdb_name,
            "\n",
            chain,
            "\n",
            sequence,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)
    elif status_code == 204:
      if response.text == "":
        return []
      print("Check why\n",
            "Only chain sequence",
            pdb_name,
            "\n",
            chain,
            "\n",
            sequence,
            "\n",
            response,
            "\n",
            response.status_code,
            "\n",
            response.text,
            "\n\n\n",
            flush=True)
  if response.status_code == 204 or response.status_code == 400:
    return []
  var_result = json.loads(response.text)

  result = []
  # Make sure that only pdb appear that we can process
  chain_score = {}
  dirty_list = []
  for i in var_result["result_set"]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      dirty_list.append(i["identifier"].split("-")[0].lower())
      chain_score[i["identifier"].split("-")[0].lower()] = i["score"]

  clean_list = np.intersect1d(np.array(dirty_list), np.array(get_all_pdb_work())).tolist()
  for i in clean_list[:can]:
    result.append([i, chain_score[i]])

  return result


def get_pdb_adn_arn_online():
  global adn_arn_online_list
  if adn_arn_online_list == []:
    adn_arn_online_list = get_pdb_adn_arn_online_aux()
  return adn_arn_online_list


def get_pdb_adn_arn_online_aux():
  search_request = {
    "query": {
      "type": "group",
      "logical_operator": "or",
      "nodes": [
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "struct_keywords.pdbx_keywords",
            "negation": False,
            "operator": "contains_phrase",
            "value": "DNA"
          }
        },
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "struct_keywords.pdbx_keywords",
            "negation": False,
            "operator": "contains_phrase",
            "value": "RNA"
          }
        },
        {
          "type": "terminal",
          "service": "text",
          "parameters": {
            "attribute": "struct_keywords.pdbx_keywords",
            "negation": False,
            "operator": "contains_phrase",
            "value": "DNA-RNA"
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "return_all_hits": True,
      "scoring_strategy": "combined",
      "sort": [
        {
          "sort_by": "score",
          "direction": "desc"
        }
      ]
    }
  }
  json_dump = json.dumps(search_request)
  url_get = 'https://search.rcsb.org/rcsbsearch/v1/query?json={0}'.format(json_dump)

  status_code = 500
  while status_code != 200:
    response = requests.get(url_get)
    status_code = response.status_code
    time.sleep(random.randint(2, 15))
    if status_code != 200 and status_code != 204:
      print(response, response.status_code, response.text, "\n\n\n", flush=True)

  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"]:
    result.append(i["identifier"].split("-")[0].lower())
  return result


def get_pdb_no_work():
  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'

  if os.path.exists(know_pdb_path):
    pd_data_frame = pd.read_csv(know_pdb_path)
    temp = pd_data_frame.values.tolist()
    result = []
    for i in temp:
      if i[1] == 0:
        result.append(i[0])
    return result
  return []


def get_pdb_adn_arn():
  result = get_pdb_adn_arn_online()
  return result


def get_ignore_pdbs():
  return get_pdb_no_work() + get_pdb_adn_arn()


def get_chains_pdb(input_file):
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
  list_result = list(dict.fromkeys(list_result))
  return list_result


def get_atoms_of_list_pdb(input_file, list_chains):
  result = {}
  for chain in list_chains:
    result[chain] = []

  input_file = os.path.abspath(input_file)

  with open(input_file) as origin_file:
    for line in origin_file:

      if line[0:4] == "ATOM":
        line_list = list(line)
        line_list[21] = "A"
        line_list = "".join(line_list)

        chain_check = line[72:76].replace(" ", "")
        if chain_check == "" or \
          chain_check is None or \
          chain_check == " " or \
          chain_check == '' or \
          chain_check == ' ':
          continue

        if chain_check in list_chains:
          result[chain_check].append(line_list)
        else:
          raise ValueError("Error in get atoms for chain, invalid chain, chain check:" + chain_check
                           + ",list_chains:" + str(list_chains)
                           + ",file:" + input_file)

  for key in result.keys():
    if result[key] == []:
      print(result)
      raise ValueError("Error in get atoms for chain, chain not exist")
  return result


def get_chain_mmcif_2_pdb(input_file):
  result = []

  input_file = os.path.abspath(input_file)

  with open(input_file) as origin_file:
    for line in origin_file:

      if line[0:4] == "ATOM":
        check_list = list(line[72:76])
        chain_id = ""
        for i in check_list:
          if i == " ":
            break
          else:
            chain_id += i

        if chain_id not in result:
          if chain_id == None or \
            chain_id == "" or \
            chain_id == " " or \
            chain_id == '' or \
            chain_id == '':
            continue
          else:
            result.append(chain_id)

  return result


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
  max_val = 2 * max(math.ceil(x_actual), math.ceil(y_actual), math.ceil(z_actual))
  max_val += 10
  return [max_val, max_val, max_val]


def move_pdb_center(pdb_path):
  input_file = os.path.abspath(pdb_path)
  atoms = []

  with open(input_file) as origin_file:
    for line in origin_file:
      if line[0:4] == "ATOM":
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atoms.append([x, y, z])

  center = np.mean(atoms, axis=0)

  all_lines = []
  with open(input_file) as origin_file2:
    for line in origin_file2:
      if line[0:4] == "ATOM":
        x = float(line[30:38]) - center[0]
        y = float(line[38:46]) - center[1]
        z = float(line[46:54]) - center[2]

        text_x = "{:8.3f}".format(x)
        text_y = "{:8.3f}".format(y)
        text_z = "{:8.3f}".format(z)

        line = change_string(30, 37, line, text_x)
        line = change_string(38, 45, line, text_y)
        line = change_string(46, 53, line, text_z)
      all_lines.append(line)

  final_text = "".join(all_lines)
  os.remove(pdb_path)
  exit_file = open(pdb_path, "w")
  exit_file.write(final_text)
  exit_file.close()


def get_pdb_chain_sequence(pdb_path, pdb_name, chain):
  input_file = os.path.abspath(pdb_path)

  list_chain_residue = []

  with open(input_file) as origin_file:
    actual_chain = ''
    for line in origin_file:
      if line[0:4] == "ATOM":
        if actual_chain == '':
          actual_chain = line[21:22]
        elif actual_chain != line[21:22]:
          actual_chain = line[21:22]

        if actual_chain == chain:
          to_add = [line[17:20].replace(" ", ""), int(line[22:26])]
          if list_chain_residue == [] or to_add != list_chain_residue[-1]:
            list_chain_residue.append(to_add)

  # Clear duplciates
  list_chain_residue = list(k for k, _ in itertools.groupby(list_chain_residue))

  list_chain_residue = sorted(list_chain_residue, key=lambda x: x[1], reverse=False)

  result = ""
  for i in list_chain_residue:
    result += i[0]

  if result == "":
    raise Exception("Chain not exist")
  return result


def get_all_pdb_name_online():
  list_servers = [["ftp.wwpdb.org", "/pub/pdb/data/structures/all/pdb/"],
                  ["ftp.rcsb.org", "/pub/pdb/data/structures/all/pdb/"],
                  ["ftp.ebi.ac.uk", "/pub/databases/pdb/data/structures/all/pdb/"],
                  ["ftp.pdbj.org", "/pub/pdb/data/structures/all/pdb/"]]

  for i in list_servers:
    try:
      ftp = FTP()
      ftp.connect(i[0])
      ftp.login()
      ftp.cwd(i[1])
      files_list = ftp.nlst()

      result = []

      for filename in files_list:
        result.append(filename[3:-7])

      return result
    except Exception as e:
      pass

  return e


def get_all_pdb_name():
  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'

  if os.path.exists(know_pdb_path):
    pd_data_frame = pd.read_csv(know_pdb_path)
    result = pd_data_frame["Name"].tolist()
    return result
  else:
    return get_all_pdb_name_online()


def get_all_pdb_work():
  from general_utils.database_utils import get_all_archive_pdb
  # all_names = get_all_pdb_name()  # 169315
  # all_names = np.setdiff1d(np.array(all_names), np.array(get_ignore_pdbs())).tolist()
  # return all_names
  return np.setdiff1d(np.array(get_all_archive_pdb()), np.array(get_ignore_pdbs())).tolist()


def get_percentage_pbs_check_file(percentage_data_set, file_checkpoint, executor, min_can_chains=0):
  from experiment.utils_general import pdb_percentage

  if not os.path.exists(file_checkpoint):
    all_names = pdb_percentage(percentage_data_set, executor, min_can_chains)  # 169315
    open_file = open(file_checkpoint, "wb")
    pickle.dump(all_names, open_file)
    open_file.close()
  else:
    open_file = open(file_checkpoint, "rb")
    all_names = pickle.load(open_file)
    open_file.close()

  return all_names


def align_pdb_file_1_in_2(pdb_file1, pdb_file2):
  from pymol import cmd
  cmd.reinitialize()

  name1 = "PDBN1"
  name2 = "PDBN2"

  cmd.load(pdb_file1, name1)
  cmd.load(pdb_file2, name2)

  try:
    result = cmd.align(name1, name2)
    finalResult = PDBsAlingResult(cmd.count_atoms(name1),
                                  cmd.count_atoms(name2),
                                  result[0],
                                  result[1],
                                  result[2],
                                  result[3],
                                  result[4],
                                  result[5],
                                  result[6])
    return finalResult
  except:
    return None
