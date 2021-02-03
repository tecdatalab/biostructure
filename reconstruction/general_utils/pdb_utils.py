import json
import math
import os
import shutil
import tempfile
from ftplib import FTP

import numpy as np
import requests


from general_utils.string_utils import change_string
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from pymol import cmd


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
  response = requests.get(url_get)
  if response.status_code == 204:
    return []
  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"][:can]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      result.append([i["identifier"].split("-")[0].lower(), i["score"]])
  return result


def get_similar_pdb_chain_structural(pdb_name, chain, can=10):
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
              "asym_id": chain.upper()
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
  response = requests.get(url_get)
  if response.status_code == 204:
    return []
  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"][:can]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      result.append([i["identifier"].split("-")[0].lower(), i["score"]])
  return result


def get_similar_pdb_chain_sequential(pdb_name, chain, can=10):
  from general_utils.download_utils import download_pdb
  path_temp = tempfile.mkdtemp()
  path_temp = os.path.abspath(path_temp)
  temp_file_path = path_temp + "/" + pdb_name + ".pdb"
  download_pdb(pdb_name, temp_file_path)
  sequence = get_pdb_chain_sequence(temp_file_path, chain)
  shutil.rmtree(path_temp)
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
                "operator": "in",
                "negation": True,
                "value": list(map(lambda x: x.upper(), get_pdb_no_work())),
                "attribute": "rcsb_entry_container_identifiers.entry_id"
              }
            },
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
            "target": "pdb_protein_sequence",
            "value": sequence
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "return_all_hits": True,
      "scoring_strategy": "sequence",
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
  response = requests.get(url_get)
  if response.status_code == 204:
    return []
  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"][:can]:
    if i["identifier"].split("-")[0].lower() != pdb_name.lower():
      result.append([i["identifier"].split("-")[0].lower(), i["score"]])
  return result

def get_pdb_adn_arn_online():
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
  response = requests.get(url_get)
  var_result = json.loads(response.text)
  result = []
  for i in var_result["result_set"]:
    result.append(i["identifier"].split("-")[0].lower())
  return result


def get_pdb_no_work():
  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_no_work.txt'
  if os.path.exists(evil_pdb_path):
    f_evil_pdb = open(evil_pdb_path)
    result = f_evil_pdb.read().splitlines()
    f_evil_pdb.close()
    return result
  return []


def get_pdb_adn_arn():
  evil_pdb_path = os.path.dirname(__file__) + '/../files/pdb_adn_arn.txt'
  result = []
  if os.path.exists(evil_pdb_path):
    f_pdb_adn_arn = open(evil_pdb_path)
    result = f_pdb_adn_arn.read().splitlines()
    f_pdb_adn_arn.close()
  else:
    result = get_pdb_adn_arn_online()
    f = open(evil_pdb_path, "w+")
    f.writelines("%s\n" % l for l in result)
    f.close()

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


def get_pdb_chain_sequence(pdb_path, chain):
  input_file = os.path.abspath(pdb_path)

  pdbparser = PDBParser()
  structure = pdbparser.get_structure(os.path.basename(input_file).split('.')[0], input_file)
  chains = {chain.id: seq1(''.join(residue.resname for residue in chain)) for chain in structure.get_chains()}

  query_chain = chains[chain]
  return query_chain


def mmCIF_to_pdb(mmCIF_path, exit_pdb_path):
  exit_pdb_path = os.path.abspath(exit_pdb_path)
  input_file = os.path.abspath(mmCIF_path)
  pdb_name = os.path.basename(input_file).split('.')[0]
  cmd.load(input_file, pdb_name)
  cmd.save(exit_pdb_path)


def get_all_pdb_name():
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
