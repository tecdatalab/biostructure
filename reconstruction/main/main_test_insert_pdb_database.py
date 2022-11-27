import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")


from general_utils.database_utils import get_chains_pdb_db, test_insert_get, clear_collection
from general_utils.mrc_uilts import get_mrc_level
from general_utils.mrc_uilts import get_mrc_to_pdb_aux

clear_collection()

pdb_in_test = "1yfq"

# test_insert_get()
#level=get_mrc_level("/home/adminuser/Downloads/1yfq.mrc")
print(get_mrc_level("/home/adminuser/Documents/1yfq_A.mrc","/home/adminuser/Downloads/1yfq.pdb"))
#print(get_mrc_to_pdb_aux(6.75,"/home/adminuser/Documents/1yfq_A.mrc"))
#chains = get_chains_pdb_db(pdb_in_test)

#print(chains)
