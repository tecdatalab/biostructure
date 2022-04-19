import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from general_utils.database_utils import clear_collection

clear_collection()


print("All database is cleaned")
