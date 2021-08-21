import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import glob, os

from tqdm import tqdm

from general_utils.terminal_utils import get_out


def downloadModelsPDB():
  get_out("rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ /work/lcastillo/allPDBS_tmp/")

  list_dirs = glob.glob("/work/lcastillo/allPDBS_tmp/*.ent.gz")

  with tqdm(total=len(list_dirs)) as pbar:
    for file in list_dirs:
      base_name = os.path.basename(file)
      base_name = base_name.split(".")[0]
      base_name = base_name[3:]

      final_path = "/work/lcastillo/allPDBS_unzip/{}.pdb".format(base_name)

      actual_file = os.path.dirname(file)
      actual_file = os.path.join(actual_file, 'pdb{0}.ent'.format(base_name))

      get_out("gunzip", "--force", file)
      get_out("mv", actual_file, final_path)

      pbar.update(1)



print("PDB")
downloadModelsPDB()
