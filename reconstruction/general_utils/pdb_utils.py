import json
import os

import requests


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
                "value": list(map(lambda x:x.lower(), get_pdb_no_work())),
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
      "pager": {
        "start": 0,
        "rows": can
      },
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
  for i in var_result["result_set"]:
    result.append([i["identifier"].split("-")[0].lower(), i["score"]])
  return result


def get_similar_pdb_chain(pdb_name, chain, can=10):
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
                "value": list(map(lambda x: x.lower(), get_pdb_no_work())),
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
              "assembly_id": chain.upper()
            },
            "operator": "relaxed_shape_match"
          }
        }
      ]
    },
    "return_type": "entry",
    "request_options": {
      "pager": {
        "start": 0,
        "rows": can
      },
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
  for i in var_result["result_set"]:
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
