import os
import datetime
import shutil
from progress.bar import Bar

folder = "temp_load_database"
minutes_delete = 60

dir_to_search = "/work/lcastillo/" + folder
all_dirs = os.listdir( dir_to_search )

bar1 = Bar('Process:', max=len(all_dirs))
print("Start")
for dirname in all_dirs:
   curpath = os.path.join(dir_to_search, dirname )
   file_modified = datetime.datetime.fromtimestamp(os.path.getmtime(curpath))
   if datetime.datetime.now() - file_modified > datetime.timedelta(minutes=minutes_delete):
      if os.path.isfile(curpath):
         os.remove(curpath)
      else:
         shutil.rmtree(curpath)
   bar1.next()
bar1.finish()
