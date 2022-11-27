import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")


from general_utils.database_utils import get_chains_pdb_db, test_insert_get, clear_collection

clear_collection()

pdb_in_test = "1yfq"

# test_insert_get()

chains = get_chains_pdb_db(pdb_in_test)



print(chains)
