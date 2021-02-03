import os

from general_utils.download_utils import download_pdb
from general_utils.pdb_utils import get_pdb_chain_sequence, get_similar_pdb_struct, get_similar_pdb_chain_structural, \
  get_similar_pdb_chain_sequential, get_chains_pdb
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains

#pdb = '5T4P'
path = "./"

pdb = '6VM4'

# download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
# chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))


#download_pdb(pdb, pdb_path)
#secuence = get_pdb_chain_sequence(pdb_path, "A")
#print(secuence)

#Todos se me parecen
# result_struct = get_similar_pdb_struct(pdb, -1)
# print(result_struct)


result_chain_struct = get_similar_pdb_chain_structural(pdb, "A")
print(result_chain_struct)
#
# for i in result_chain_struct:
#   pdb = i[0]
#
#   download_pdb(pdb, '{0}/{1}.pdb'.format(path, pdb))
#   chains = get_chains_pdb('{0}/{1}.pdb'.format(path, pdb))
#   pdb_to_mrc_chains(False, False, 5.0, '{0}/{1}.pdb'.format(path, pdb), path, chains,
#                     len(chains))
#   os.remove('{0}/{1}.pdb'.format(path, pdb))


# result_chain_sequence = get_similar_pdb_chain_sequential(pdb, "A")
# print(result_chain_sequence)
