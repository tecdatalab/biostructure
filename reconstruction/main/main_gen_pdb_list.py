import sys
import pathlib

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
import os
import shutil
import tempfile
import urllib.request
import pandas as pd
import numpy as np
import random
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from general_utils.download_utils import download_pdb
from general_utils.pdb_utils import get_chains_pdb
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_mrc.generate import get_mrc_synthetic_segments_pdb


def test_pdb(pdb_name):
  pdb_name = pdb_name.lower()
  print("\n\n\n Enter:" + pdb_name + "\n\n\n", flush=True)
  if os.path.exists("/work/lcastillo"):
    path_dir = os.path.abspath(tempfile.mkdtemp(dir="/work/lcastillo/temp"))
  else:
    path_dir = os.path.abspath(tempfile.mkdtemp())

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(path_dir, pdb_name))
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_dir, pdb_name))

  pdb_to_mrc_chains(True, False, 5.0, '{0}/{1}.pdb'.format(path_dir, pdb_name), path_dir, chains, len(chains))
  local_path = path_dir + "/" + pdb_name

  _segments, _original_structure = get_mrc_synthetic_segments_pdb('{0}/{1}.mrc'.format(local_path, pdb_name),
                                                                  local_path, calculate_Z3D=True)

  shutil.rmtree(path_dir)

  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'
  write_in_file(know_pdb_path, ["Name", "OK"], [[pdb_name, 1]])
  print("\n\n\n Finish:" + pdb_name + "\n\n\n", flush=True)


def gen_update_pdb_list():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:
      if os.path.exists("/work/lcastillo"):
        if os.path.exists("/work/lcastillo/temp"):
          shutil.rmtree("/work/lcastillo/temp")
        os.mkdir("/work/lcastillo/temp")

      path_dir = os.path.abspath(tempfile.mkdtemp())
      urllib.request.urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/author.idx", path_dir + "/data.txt")
      with open(path_dir + "/data.txt") as f:
        content = f.readlines()
      real_pdb_name = []
      for i in content:
        if i.find(" ;")!=-1:
          split_data = i.split(" ;")
          if len(split_data[0])==4:
            real_pdb_name.append(split_data[0].lower())
      shutil.rmtree(path_dir)


      #real_pdb_name = ["5wob"]
      know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'

      if os.path.exists(know_pdb_path):
        pd_data_frame = pd.read_csv(know_pdb_path)
        actual_pdb_list = pd_data_frame["Name"].tolist()
      else:
        actual_pdb_list = []

      real_pdb_name = np.setdiff1d(real_pdb_name, actual_pdb_list).tolist()
      #real_pdb_name = ["5uyk"]
      print(len(real_pdb_name))
      random.shuffle(real_pdb_name)

      #real_pdb_name = ["1fnt"]

      parallel_jobs = []
      for pdb_name in real_pdb_name:
        parallel_jobs.append([pdb_name, executor.submit(test_pdb, pdb_name)])
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          print(f[0])
          print(e)
          write_in_file(know_pdb_path, ["Name", "OK"], [[f[0], 0]])


if __name__ == '__main__':
  gen_update_pdb_list()
