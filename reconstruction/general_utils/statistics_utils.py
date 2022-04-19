import os
import pandas as pd


def get_dicc_pdbs_can_chains_file():
  result = {}
  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_can_chains.csv'
  df = pd.read_csv(know_pdb_path)
  data = df.values.tolist()

  for pdb_data in data:
    list_pdbs = result.get(pdb_data[1], [])
    list_pdbs.append(pdb_data[0])
    result[pdb_data[1]] = list_pdbs

  return result
