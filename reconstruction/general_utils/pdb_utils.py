import json
import requests


def get_similar_pdb_struct(pdb_name, can=10):
  search_request = {
    "query": {
      "type": "terminal",
      "service": "structure",
      "parameters": {
        "value": {
          "entry_id": pdb_name.upper(),
          "assembly_id": "1"
        },
        "operator": "relaxed_shape_match"
      }
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
      "type": "terminal",
      "service": "structure",
      "parameters": {
        "value": {
          "entry_id": pdb_name.upper(),
          "asym_id": chain.upper(),
        },
        "operator": "relaxed_shape_match"
      }
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
