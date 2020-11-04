import os
import urllib.request
import progressbar
import shutil
from general_utils.terminal_utils import get_out
import tempfile


class MyProgressBar:
  def __init__(self):
    self.pbar = None

  def __call__(self, block_num, block_size, total_size):
    if not self.pbar:
      self.pbar = progressbar.ProgressBar(maxval=total_size)
      self.pbar.start()

    downloaded = block_num * block_size
    if downloaded < total_size:
      self.pbar.update(downloaded)
    else:
      self.pbar.finish()


def download_emd(id_code, exit_path, create_progress_bar=False):
  # path = './temp_emd_download'
  path = tempfile.mkdtemp()
  path = os.path.abspath(path)
  exit_path = os.path.abspath(exit_path)
  file_url = 'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{0}/map/emd_{0}.map.gz'.format(id_code)

  if os.path.exists(path):
    shutil.rmtree(path)
  os.mkdir(path)

  if create_progress_bar:
    urllib.request.urlretrieve(file_url, path + '/emd_{0}.map.gz'.format(id_code), MyProgressBar())
  else:
    urllib.request.urlretrieve(file_url, path + '/emd_{0}.map.gz'.format(id_code))

  get_out("gunzip", path + '/emd_{0}.map.gz'.format(id_code))
  get_out("rm", path + '/emd_{0}.map.gz'.format(id_code))
  get_out("mv", path + '/emd_{0}.map'.format(id_code), exit_path)

  shutil.rmtree(path)


def download_pdb(id_code, exit_path, create_progress_bar=False):
  # path = './temp_emd_download'
  path = tempfile.mkdtemp()
  path = os.path.abspath(path)
  exit_path = os.path.abspath(exit_path)
  file_url = 'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/pdb/pdb{0}.ent.gz'.format(id_code)

  if os.path.exists(path):
    shutil.rmtree(path)
  os.mkdir(path)

  if create_progress_bar:
    urllib.request.urlretrieve(file_url, path + '/pdb{0}.ent.gz'.format(id_code), MyProgressBar())
  else:
    urllib.request.urlretrieve(file_url, path + '/pdb{0}.ent.gz'.format(id_code))

  get_out("gunzip", path + '/pdb{0}.ent.gz'.format(id_code))
  get_out("mv", path + '/pdb{0}.ent'.format(id_code), exit_path)

  shutil.rmtree(path)


def get_all_pdb_name():
  from ftplib import FTP
  ftp = FTP()
  ftp.connect("ftp.ebi.ac.uk")
  ftp.login()
  ftp.cwd("/pub/databases/pdb/data/structures/all/pdb/")
  files_list = ftp.nlst()

  result = []

  for filename in files_list:
    result.append(filename[3:-7])

  return result


def get_all_emd_name():
  from ftplib import FTP
  ftp = FTP()
  ftp.connect("ftp.wwpdb.org")
  ftp.login()
  ftp.cwd("/pub/emdb/structures/")
  files_list = ftp.nlst()

  result = []

  for filename in files_list:
    result.append(filename[4:])

  return result
