import os
import urllib.request
import progressbar
import shutil

from general_utils.temp_utils import gen_dir, free_dir
from general_utils.terminal_utils import get_out
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
    er = None
    for url_format_string in url_format_principal_list:
      try:
        download_biomolecular_zip_file(id_code, exit_path, url_format_string, 'emd_{0}.map.gz', 'emd_{0}.map',
                                       create_progress_bar)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_emd_xml(id_code, exit_path, create_progress_bar=False):
  can_try = 60
  url_format_principal_list = ['ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-{0}/header/emd-{0}.xml',
                               'ftp://ftp.rcsb.org/pub/emdb/structures/EMD-{0}/header/emd-{0}.xml',
                               'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{0}/header/emd-{0}.xml',
                               'ftp://ftp.pdbj.org/pub/emdb/structures/EMD-{0}/header/emd-{0}.xml']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    er = None
    for url_format_string in url_format_principal_list:
      try:
        download_biomolecular_file(id_code, exit_path, url_format_string, create_progress_bar)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_pdb(id_code, exit_path, create_progress_bar=False):
  id_code = id_code.lower()
  ok_flag = False
  er = None
  try:
    download_pdb_in_pdb(id_code, exit_path, create_progress_bar)
    ok_flag = True
  except Exception as e:
    er = e
    ok_flag = False

  if not ok_flag:
    try:
      download_pdb_in_http(id_code, exit_path, create_progress_bar)
      ok_flag = True
    except Exception as e:
      er = e
      ok_flag = False

  if not ok_flag:
    try:
      download_pdb_in_mmCIF(id_code, exit_path, create_progress_bar)
      ok_flag = True
    except Exception as e:
      er = e
      ok_flag = False

  if not ok_flag:
    try:
      download_pdb_in_mmCIF_http(id_code, exit_path, create_progress_bar)
      ok_flag = True
    except Exception as e:
      er = e
      ok_flag = False

  if not ok_flag:
    raise er
  else:
    return True


def download_pdb_in_http(id_code, exit_path, create_progress_bar=False):
  can_try = 30
  url_format_principal_list = ['https://files.rcsb.org/download/{0}.pdb']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    er = None
    for url_format_string in url_format_principal_list:
      try:
        download_biomolecular_file(id_code.upper(), exit_path, url_format_string, create_progress_bar)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_pdb_in_mmCIF_http(id_code, exit_path, create_progress_bar=False):
  from general_utils.pdb_utils import mmCIF_to_pdb
  er = None
  can_try = 30
  url_format_principal_list = ['https://files.rcsb.org/download/{0}.cif']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    for url_format_string in url_format_principal_list:
      try:
        path_temp = gen_dir()
        temp_file_path = path_temp + "/" + os.path.basename(exit_path).split(".")[0] + ".cif"
        download_biomolecular_file(id_code.upper(), temp_file_path, url_format_string, create_progress_bar)
        mmCIF_to_pdb(temp_file_path, exit_path)
        free_dir(path_temp)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_pdb_in_pdb(id_code, exit_path, create_progress_bar=False):
  can_try = 30
  url_format_principal_list = ['ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.rcsb.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/pdb/pdb{0}.ent.gz',
                               'ftp://ftp.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    er = None
    for url_format_string in url_format_principal_list:
      try:
        download_biomolecular_zip_file(id_code, exit_path, url_format_string, 'pdb{0}.ent.gz', 'pdb{0}.ent',
                                       create_progress_bar)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_pdb_in_mmCIF(id_code, exit_path, create_progress_bar=False):
  from general_utils.pdb_utils import mmCIF_to_pdb
  er = None
  can_try = 30
  url_format_principal_list = ['ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/{0}.cif.gz',
                               'ftp://ftp.rcsb.org/pub/pdb/data/structures/all/mmCIF/{0}.cif.gz',
                               'ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/mmCIF/{0}.cif.gz',
                               'ftp://ftp.pdbj.org/pub/pdb/data/structures/all/mmCIF/{0}.cif.gz']
  while True:
    random.shuffle(url_format_principal_list)
    flag = False
    for url_format_string in url_format_principal_list:
      try:
        path_temp = gen_dir()
        temp_file_path = path_temp + "/" + os.path.basename(exit_path).split(".")[0] + ".cif"
        download_biomolecular_zip_file(id_code, temp_file_path, url_format_string, '{0}.cif.gz', '{0}.cif',
                                       create_progress_bar)
        mmCIF_to_pdb(temp_file_path, exit_path)
        free_dir(path_temp)
        flag = True
        break
      except Exception as e:
        er = e
        permit_errors_codes = ["421", "104"]
        flag_error = False
        for i in permit_errors_codes:
          if str(e).find(i) != -1:
            flag_error = True
            break
        if not flag_error:
          raise e

    if not flag:
      can_try -= 1
      if can_try < 0:
        raise er
      else:
        time.sleep(15)
    else:
      break


def download_biomolecular_zip_file(id_code, exit_path, url_format_string, file_zip, file_unzip,
                                   create_progress_bar=False):
  # path = './temp_emd_download'
  path = gen_dir()
  path = os.path.abspath(path)
  exit_path = os.path.abspath(exit_path)
  file_url = url_format_string.format(id_code)

  if os.path.exists(path):
    shutil.rmtree(path)
  os.mkdir(path)

  ###-------------
  try:
    with open(path + "/" + file_zip.format(id_code), 'w') as fp:
      pass
  except:
    pass
  ###------------

  if create_progress_bar:
    urllib.request.urlretrieve(file_url, path + "/" + file_zip.format(id_code), MyProgressBar())
  else:
    urllib.request.urlretrieve(file_url, path + "/" + file_zip.format(id_code))

  get_out("gunzip", path + "/" + file_zip.format(id_code))
  get_out("mv", path + "/" + file_unzip.format(id_code), exit_path)

  free_dir(path)


def download_biomolecular_file(id_code, exit_path, url_format_string, create_progress_bar=False):
  # path = './temp_emd_download'
  exit_path = os.path.abspath(exit_path)
  file_url = url_format_string.format(id_code)

  ###-----------------
  try:
    with open(exit_path, 'w') as fp:
      pass
  except:
    pass
  ###-----------------

  if create_progress_bar:
    urllib.request.urlretrieve(file_url, exit_path, MyProgressBar())
  else:
    urllib.request.urlretrieve(file_url, exit_path)
