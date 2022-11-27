import csv  # imports module csv
import os
import sys
import pathlib
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

from general_utils.temp_utils import gen_dir, free_dir
from general_utils.download_utils import download_pdb_in_pdb


def exist_pdb_file(pdb_id, folder_work):
  exit_file = "{}/{}.pdb".format(folder_work, pdb_id)
  download_pdb_in_pdb(pdb_id, exit_file)
  if os.path.getsize(exit_file) == 0:
    return False
  os.remove(exit_file)
  return True

def main():

  # Parale
  comm = MPI.COMM_WORLD
  size = comm.Get_size()

  with MPICommExecutor(comm, root=0, worker_size=size) as executor:
    if executor is not None:

      know_pdb_path = os.path.dirname(__file__) + '/../files/pdb_can_chains.csv'
      output = os.path.dirname(__file__) + '/../files/pdb_can_chains_with_pdb_flag.csv'
      delim = ","  # set your own delimiter

      if not os.path.exists(output):
        file_h = open(output, "w+")
        file_h.close()

      source1 = csv.reader(open(know_pdb_path, "r"), delimiter=delim)
      source2 = csv.reader(open(know_pdb_path, "r"), delimiter=delim)
      output_source1 = csv.reader(open(output, "r"), delimiter=delim)

      row0 = next(source1)
      row0.append("ExistPDBFile")
      # open csv readers

      source_dict = {}
      for row in output_source1:
        source_dict[row[0]] = True

      total_process = sum(1 for line in source2)
      actual_process = 0

      folder_work = gen_dir()

      #Create Works
      parallel_jobs = []
      all_keys = source_dict.keys()
      for row in source1:
        if not row[0] in all_keys:
          parallel_jobs.append([row[0], row, executor.submit(exist_pdb_file, row[0], folder_work)])

      # write new changed rows
      with open(output, "a+") as fout:
        csvwriter = csv.writer(fout, delimiter=delim)

        for f in parallel_jobs:
          try:
            has_pdb_file = f[2].result()
            row = f[1]
            row.append(has_pdb_file)
            csvwriter.writerow(row)
          except Exception as e:
            print(f[0], e, flush=True)

          actual_process+=1
          print( "{}%".format((actual_process/total_process)*100), flush=True)

      free_dir(folder_work)

main()
