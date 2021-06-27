import os
import datetime
import shutil
from tqdm import tqdm
import glob

listing = glob.glob("/work/lcastillo/data_experiment_1_c*")
for filename in tqdm(listing):
  if os.path.isdir(filename):
    shutil.rmtree(filename)
