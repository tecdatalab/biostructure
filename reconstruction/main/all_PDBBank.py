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
  from constans.directory_constans import CLUSTER_ZIPPED_PDBS_FOLDER, CLUSTER_UNZIPPED_PDBS_FOLDER
  zipped_folder = CLUSTER_ZIPPED_PDBS_FOLDER
  unzipped_folder = CLUSTER_UNZIPPED_PDBS_FOLDER

  # exist or not.
  if not os.path.isdir(zipped_folder):
    os.makedirs(zipped_folder)

  if not os.path.isdir(unzipped_folder):
    os.makedirs(unzipped_folder)

  execute_command("rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ {}/".format(zipped_folder))

  for root, dirs, files in tqdm(os.walk("{}/".format(zipped_folder))):
    for file in files:
      if file.endswith('.ent.gz'):
        file = os.path.join(root, file)

        base_name = os.path.basename(file)
        base_name = base_name.split(".")[0]
        base_name = base_name[3:]

        final_path = "{}/{}.pdb".format(unzipped_folder, base_name)

        if not (os.path.exists(final_path)):

          actual_file = os.path.dirname(file)
          actual_file = os.path.join(actual_file, 'pdb{0}.ent'.format(base_name))

          commad = "gunzip -c --force " + file + " > " + final_path
          execute_command(commad)


print("PDB")
downloadModelsPDB()
