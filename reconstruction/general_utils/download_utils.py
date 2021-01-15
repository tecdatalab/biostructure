import os
import urllib.request
from urllib.error import URLError
import progressbar
import shutil
from general_utils.terminal_utils import get_out
import tempfile
import time
import random


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
  can_try = 60
  url_format_principal_list = ['ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz',
                               'ftp://ftp.rcsb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz',
                               'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{0}/map/emd_{0}.map.gz',
                               'ftp://ftp.pdbj.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    for url_format_string in url_format_principal_list:
      try:
        download_emd_aux(id_code, exit_path, url_format_string, create_progress_bar)
        flag = True
        break
      except Exception as e:
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          return e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise e
      else:
        time.sleep(15)
    else:
      break


def download_emd_aux(id_code, exit_path, url_format_string, create_progress_bar=False):
  # path = './temp_emd_download'
  path = tempfile.mkdtemp()
  path = os.path.abspath(path)
  exit_path = os.path.abspath(exit_path)
  file_url = url_format_string.format(id_code)

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
  can_try = 60
  url_format_principal_list = ['ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    for url_format_string in url_format_principal_list:
      try:
        download_pdb_aux(id_code, exit_path, url_format_string, create_progress_bar)
        flag = True
        break
      except Exception as e:
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          return e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise e
      else:
        time.sleep(15)
    else:
      break


def download_pdb_aux(id_code, exit_path, url_format_string, create_progress_bar=False):
  # path = './temp_emd_download'
  path = tempfile.mkdtemp()
  path = os.path.abspath(path)
  exit_path = os.path.abspath(exit_path)
  file_url = url_format_string.format(id_code)

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
  list_servers = [["ftp.wwpdb.org", "/pub/pdb/data/structures/all/pdb/"],
                  ["ftp.rcsb.org", "/pub/pdb/data/structures/all/pdb/"],
                  ["ftp.ebi.ac.uk", "/pub/databases/pdb/data/structures/all/pdb/"],
                  ["ftp.pdbj.org", "/pub/pdb/data/structures/all/pdb/"]]

  for i in list_servers:
    try:
      ftp = FTP()
      ftp.connect(i[0])
      ftp.login()
      ftp.cwd(i[1])
      files_list = ftp.nlst()

      result = []

      for filename in files_list:
        result.append(filename[3:-7])

      return result
    except Exception as e:
      pass

  return e


def get_all_emd_name():
  from ftplib import FTP
  list_servers = [["ftp.wwpdb.org", "/pub/emdb/structures/"],
                  ["ftp.rcsb.org", "/pub/emdb/structures/"],
                  ["ftp.ebi.ac.uk", "/pub/databases/emdb/structures/"],
                  ["ftp.pdbj.org", "/pub/emdb/structures/"]]

  for i in list_servers:
    try:
      ftp = FTP()
      ftp.connect(i[0])
      ftp.login()
      ftp.cwd(i[1])
      files_list = ftp.nlst()

      result = []

      for filename in files_list:
        result.append(filename[4:])

      return result
    except Exception as e:
      pass
  return e
