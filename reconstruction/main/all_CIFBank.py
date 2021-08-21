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

def downloadModelsCIF():
  execute_command("rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ /work/lcastillo/allCIFS_tmp/")

  for root, dirs, files in tqdm(os.walk("/work/lcastillo/allCIFS_tmp/")):
    for file in files:
      if file.endswith('.cif.gz'):
        file = os.path.join(root, file)

        base_name = os.path.basename(file)
        base_name = base_name.split(".")[0]

        final_path = "/work/lcastillo/allCIFS_unzip/{}.cif".format(base_name)

        if not (os.path.exists(final_path)):

          actual_file = os.path.dirname(file)
          actual_file = os.path.join(actual_file, 'pdb{0}.cif'.format(base_name))

          get_out("gunzip", "--keep", "--force", file)
          get_out("mv", actual_file, final_path)



print("CIF")
downloadModelsCIF()


