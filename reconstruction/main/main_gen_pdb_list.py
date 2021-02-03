import os
import shutil
import tempfile
import urllib.request
import pandas as pd
import numpy as np
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

from csv_modules.csv_writer import write_in_file
from general_utils.download_utils import download_pdb
from general_utils.pdb_utils import get_chains_pdb
from pdb_to_mrc.pdb_2_mrc import pdb_to_mrc_chains
from process_mrc.generate import get_mrc_synthetic_segments_pdb


def test_pdb(pdb_name):
  pdb_name = pdb_name.lower()
  path_dir = os.path.abspath(tempfile.mkdtemp())

  download_pdb(pdb_name, '{0}/{1}.pdb'.format(path_dir, pdb_name))
  chains = get_chains_pdb('{0}/{1}.pdb'.format(path_dir, pdb_name))

  pdb_to_mrc_chains(True, False, 5.0, '{0}/{1}.pdb'.format(path_dir, pdb_name), path_dir, chains, len(chains))
  local_path = path_dir+"/"+pdb_name

  _segments, _original_structure = get_mrc_synthetic_segments_pdb('{0}/{1}.mrc'.format(local_path, pdb_name),
                                                                       local_path, calculate_Z3D=True)

  shutil.rmtree(path_dir)

  know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'
  write_in_file(know_pdb_path, ["Name", "OK"], [[pdb_name, 1]])


def gen_update_pdb_list():
  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      f = tempfile.TemporaryDirectory()
      path_dir = os.path.abspath(tempfile.mkdtemp())
      urllib.request.urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/author.idx", path_dir + "/data.txt")
      with open(path_dir + "/data.txt") as f:
        content = f.readlines()
      real_pdb_name = []
      for i in content:
        if i.find(";")!=-1:
          split_data = i.split(";")
          if len(split_data[0])==4:
            real_pdb_name.append(split_data[0].lower())
      shutil.rmtree(path_dir)

      know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_list.csv'

      if os.path.exists(know_pdb_path):
        pd_data_frame = pd.read_csv(know_pdb_path, index_col="Name")
        actual_pdb_list = pd_data_frame.values.tolist()
      else:
        actual_pdb_list = []

      real_pdb_name = np.setdiff1d(real_pdb_name, actual_pdb_list).tolist()

      parallel_jobs = []
      for pdb_name in real_pdb_name:
        parallel_jobs.append([pdb_name, executor.submit(test_pdb, pdb_name)])
        # do_parallel_test_a_aux(path, pdb_name, result_cvs_file, resolution, range_incompleteness, can_try_experiments)
      for f in parallel_jobs:
        try:
          f[1].result()
        except Exception as e:
          print(f[0])
          print(e)

if __name__ == '__main__':
  gen_update_pdb_list()
