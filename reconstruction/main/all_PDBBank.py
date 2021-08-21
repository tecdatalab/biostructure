import subprocess
import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

import glob, os

from tqdm import tqdm

from general_utils.terminal_utils import get_out

def execute_command(cmd):
  try:
    subprocess.check_call([cmd], shell=True)
  except subprocess.CalledProcessError:
    raise RuntimeError('Command "%s" does not work' % cmd)
  except OSError:
    raise Exception('Command "%s" does not exist' % cmd)


def downloadModelsPDB():
  execute_command("rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ /work/lcastillo/allPDBS_tmp/")

  for root, dirs, files in tqdm(os.walk("/work/lcastillo/allPDBS_tmp/")):
    for file in files:
      if file.endswith('.ent.gz'):
        file = os.path.join(root, file)

        base_name = os.path.basename(file)
        base_name = base_name.split(".")[0]
        base_name = base_name[3:]

        final_path = "/work/lcastillo/allPDBS_unzip/{}.pdb".format(base_name)

        if not (os.path.exists(final_path)):

          actual_file = os.path.dirname(file)
          actual_file = os.path.join(actual_file, 'pdb{0}.ent'.format(base_name))

          get_out("gunzip", "-k", "--force", file)
          get_out("mv", actual_file, final_path)




print("PDB")
downloadModelsPDB()



