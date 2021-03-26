import os
from ftplib import FTP
import xml.etree.ElementTree as ET

from general_utils.download_utils import download_emd_xml
from general_utils.temp_utils import gen_dir, free_dir


def get_all_emd_name():
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


def get_associated_pdb(id_emd):
  path = gen_dir()
  path = os.path.abspath(path)
  result = []
  download_emd_xml(id_emd, os.path.join(path, id_emd+".xml"))

  tree = ET.parse(os.path.join(path, id_emd+".xml"))
  root = tree.getroot()
  for fittedPDB in root.iter('fittedPDBEntryId'):
    result.append(fittedPDB.text)

  free_dir(path)
  return result


