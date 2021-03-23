import os
import datetime
import shutil
from progress.bar import Bar

folder = "temp_load_database"
minutes_delete = 60

dir_to_search = "/work/lcastillo/" + folder

bar1 = Bar('Process:', max=len(list(os.walk(dir_to_search))))
print("Start")
for dirpath, dirnames, filenames in os.walk(dir_to_search):
   for file in filenames:
      curpath = os.path.join(dirpath, file)
      file_modified = datetime.datetime.fromtimestamp(os.path.getmtime(curpath))
      if datetime.datetime.now() - file_modified > datetime.timedelta(minutes=minutes_delete):
          if os.path.isfile(curpath):
            os.remove(curpath)
          else:
            shutil.rmtree(curpath)
   bar1.next()
bar1.finish()
