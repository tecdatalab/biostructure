from general_utils.download_utils import download_pdb
from general_utils.pdb_utils import get_pdb_chain_sequence

pdb = '6n2z'
pdb_path = './{0}.pdb'.format(pdb)

download_pdb(pdb, pdb_path)
secuence = get_pdb_chain_sequence(pdb_path, "A")
print(secuence)
